function processedFaaFile = processProteinFastaFile(faaFile, geneTable, headerCol, varargin)
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
%
% Name-Value Arguments
% --------------------
% outputDir : char
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

    p = parseRAVENargs(varargin, {'outputDir', pwd});
    outputDir = char(p.outputDir);

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


    % Read the FASTA, rename headers, write the result
    fastaStruct = readFasta(faaFile);
    nMatched    = 0;
    for i = 1:numel(fastaStruct)
        % The protein accession is the first whitespace-delimited token of
        % the header (readFasta strips the leading '>')
        tokens     = strsplit(fastaStruct(i).Header);
        proteinAcc = tokens{1};
        if isKey(proteinMap, proteinAcc)
            fastaStruct(i).Header = headerValues{proteinMap(proteinAcc)};
            nMatched = nMatched + 1;
        end
        % Sequences with no match keep their original header
    end
    writeFasta(processedFaaFile, fastaStruct);

    fprintf('Processed FASTA written to: %s\n', processedFaaFile);
    fprintf('  Total sequences  : %d\n', numel(fastaStruct));
    fprintf('  Renamed          : %d\n', nMatched);
    fprintf('  Kept (no match)  : %d\n', numel(fastaStruct) - nMatched);
end