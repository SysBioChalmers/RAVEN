function processedFaaFile = processProteinFastaFile(faaFile, geneTable, headerCol, outputDir)
% processProteinFastaFile  Rename protein FASTA headers using a gene mapping table.
%
% Reads a protein FASTA file (as downloaded by downloadGenomeData) and
% replaces each sequence header with the value from the specified geneTable
% column, matched via the GenBank_protein accession present in the original
% FASTA header.  Sequences whose accession is not found in geneTable are
% kept with their original header unchanged.
%
% Parameters
% ----------
% faaFile : char
%     Path to the protein FASTA file (.faa) to process.
% geneTable : char | table
%     Either a path to a tab-delimited .tsv file produced by getGeneData,
%     or a MATLAB table variable.  Must contain at least the columns
%     'GenBank_protein' and the column named by headerCol.
% headerCol : char
%     Name of the geneTable column whose values will replace each FASTA
%     header (e.g. 'locus_tag', 'gene_name', 'GenBank_protein').
% outputDir : char, optional
%     Directory where the processed FASTA file is saved.  The output file
%     name is the original base name with '_processed' appended before the
%     extension (default: current working directory).
%
% Returns
% -------
% processedFaaFile : char
%     Full path to the written processed FASTA file.
%
% Examples
% --------
%     % Rename headers using locus_tag
%     [~, faa]  = downloadGenomeData('GCF_000002595.2');
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv');
%     outFile   = processProteinFastaFile(faa, geneTable, 'locus_tag');
%
%     % Load both inputs from saved files, use gene_name as header
%     outFile = processProteinFastaFile('GCF_000002595.2_protein.faa', 'chlamy.tsv', 'gene_name');
%
%     % Save the processed file to a specific directory
%     outFile = processProteinFastaFile(faa, 'chlamy.tsv', 'locus_tag', 'results/');
%
% Notes
% -----
% GenBank_protein must contain accessions matching the first token of each
% FASTA header (after '>').  For prokaryotes this is the protein_id
% (e.g. WP_012345678.1); for eukaryotes it is the GenBank accession stored
% in the GFF3 annotation by getGeneData.
	
    if nargin < 1 || isempty(faaFile)
        error('processProteinFastaFile:missingInput', ...
            'faaFile is required. Provide a path to a protein FASTA (.faa) file.');
    end
    faaFile = char(faaFile);
    if ~isfile(faaFile)
        error('processProteinFastaFile:fileNotFound', ...
            'FASTA file not found: ''%s''.', faaFile);
    end

    if nargin < 2 || isempty(geneTable)
        error('processProteinFastaFile:missingInput', ...
            'geneTable is required. Provide a table variable or a path to a .tsv file.');
    end

    if nargin < 3 || isempty(headerCol)
        error('processProteinFastaFile:missingInput', ...
            'headerCol is required. Specify the geneTable column to use as the new FASTA header (e.g. ''locus_tag'').');
    end
    headerCol = char(headerCol);

    if nargin < 4 || isempty(outputDir)
        outputDir = pwd;
    else
        outputDir = char(outputDir);
    end

    % Load gene table
    if ischar(geneTable) || isstring(geneTable)
        tsvPath = char(geneTable);
        if ~isfile(tsvPath)
            error('processProteinFastaFile:fileNotFound', ...
                'Gene table file not found: ''%s''.', tsvPath);
        end
        geneTable = readtable(tsvPath, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'char');
    end

    if ~istable(geneTable)
        error('processProteinFastaFile:badInput', ...
            'geneTable must be a MATLAB table or a path to a .tsv file.');
    end

    for col = {'GenBank_protein', headerCol}
        if ~ismember(col{1}, geneTable.Properties.VariableNames)
            error('processProteinFastaFile:missingColumn', ...
                'geneTable is missing required column: ''%s''.', col{1});
        end
    end

    % Build protein-accession → header-value lookup map
    proteinIDs = geneTable.GenBank_protein;
    if ~iscell(proteinIDs), proteinIDs = cellstr(proteinIDs); end

    headerValues = geneTable.(headerCol);
    if ~iscell(headerValues), headerValues = cellstr(headerValues); end

    % Keep only rows with a non-empty protein accession
    validMask  = ~cellfun(@isempty, proteinIDs);
    proteinMap = containers.Map(proteinIDs(validMask), num2cell(find(validMask)));


    % Prepare output path
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    [~, baseName, ext] = fileparts(faaFile);
    processedFaaFile   = fullfile(outputDir, [baseName '_processed' ext]);


    % Stream FASTA, rename headers, write output
    nTotal     = 0;
    nMatched   = 0;
    nUnmatched = 0;

    fidIn  = fopen(faaFile,          'r');
    fidOut = fopen(processedFaaFile, 'w');

    currentHeader = '';
    seqBuffer     = {};

    while ~feof(fidIn)
        line = fgetl(fidIn);
        if ~ischar(line), break, end

        if startsWith(line, '>')
            % Flush the previous record before processing the new header
            if ~isempty(currentHeader)
                writeRecord(fidOut, currentHeader, seqBuffer);
            end
            seqBuffer = {};
            nTotal    = nTotal + 1;

            % Extract the protein accession: first whitespace-delimited token
            tokens     = strsplit(strtrim(line(2:end)));
            proteinAcc = tokens{1};

            if isKey(proteinMap, proteinAcc)
                idx           = proteinMap(proteinAcc);
                currentHeader = sprintf('>%s', headerValues{idx});
                nMatched      = nMatched + 1;
            else
                currentHeader = line;   % Keep original header unchanged
                nUnmatched    = nUnmatched + 1;
            end

        else
            seqBuffer{end+1} = line; %#ok<AGROW>
        end
    end

    % Flush the last record
    if ~isempty(currentHeader)
        writeRecord(fidOut, currentHeader, seqBuffer);
    end

    fclose(fidIn);
    fclose(fidOut);

    fprintf('Processed FASTA written to: %s\n', processedFaaFile);
    fprintf('  Total sequences  : %d\n', nTotal);
    fprintf('  Renamed          : %d\n', nMatched);
    fprintf('  Kept (no match)  : %d\n', nUnmatched);
end

%--------------------------------------------------------------------------
% Helper function

function writeRecord(fid, header, lines)
% Write one FASTA record (header + sequence lines) to an open file.
    fprintf(fid, '%s\n', header);
    for k = 1:numel(lines)
        fprintf(fid, '%s\n', lines{k});
    end
end