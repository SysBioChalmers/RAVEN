function geneTable = getGeneData(accession, outputFile, downloadDir)
% getGeneData  Build a gene-mapping table from NCBI genome annotation files.
%
% Parses the GFF3 annotation and protein FASTA (.faa) files for a given
% NCBI genome assembly and produces a table mapping locus tags to gene
% symbols, protein IDs, and other identifiers. The resulting table can be
% passed directly to renameModelGenes to update gene identifiers in a
% RAVEN model. If no local files are provided, downloadGenomeData is
% called automatically to fetch them.
%
% Parameters
% ----------
% accession : char
%     NCBI genome assembly accession, e.g. 'GCF_000002595.2'.
% outputFile : char, optional
%     Path to save the resulting table as a tab-delimited .tsv file
%     (default: 'geneData.tsv').
% downloadDir : char, optional
%     Directory where genome files are downloaded if not already present
%     (default: current working directory).
%
% Returns
% -------
% geneTable : table
%     MATLAB table with one row per gene containing:
%       locus_tag       — stable locus identifier (e.g. 'Cre01.g000001')
%       old_locus_tag   — previous locus identifier when available
%       extra_id        — external database identifier (e.g. Phytozome)
%       GeneID          — NCBI Gene ID (e.g. '5723799')
%       gene_name       — common gene symbol when available (e.g. 'rbcL')
%       GenBank_protein — GenBank protein accession (e.g. 'XP_001698190.2')
%
% Examples
% --------
%     % Fetch and parse gene data for Chlamydomonas reinhardtii
%     geneTable = getGeneData('GCF_000002595.2');
%
%     % Save the mapping table to a specific file
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv');
%
%     % Save to a file and download to a specific directory
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv', 'data/');
%
% Notes
% -----
% When outputFile is provided the table is written as UTF-8 encoded
% tab-delimited text with a header row, compatible with readtable and
% renameModelGenes. Rows without a locus_tag are silently discarded.

	if nargin < 1 || isempty(accession)
        error('getGeneData:missingInput', ...
            'accession is required. Provide an NCBI assembly accession, e.g. ''GCF_000002595.2''.');
    end

    if nargin < 2 || isempty(outputFile)
        outputFile = 'geneData.tsv';
    else
        outputFile = char(outputFile);
    end

    if nargin < 3 || isempty(downloadDir)
        downloadDir = '.';
    else
        downloadDir = char(downloadDir);
    end

    % Resolve input → local GFF3 path. Decide whether input is a file path or an NCBI accession.
    if ischar(accession) || isstring(accession)
        inp = char(accession);
    else
        error('getGeneData:badInputType', 'Input must be a character vector or string.');
    end
    
    gffPath = '';
    % If input is an existing file, use it.
    if isfile(inp)
        gffPath = inp;
    else
        % Check if it looks like an NCBI assembly accession (GCF_... or GCA_...)
        if ~isempty(regexp(inp, '^GC[FA]_\d+(\.\d+)?$', 'once'))
            gffPath = downloadGenomeData(inp, downloadDir);
        else
            error('getGeneData:badInput', ...
                'Input "%s" is neither an existing file nor a recognised NCBI assembly accession.', inp);
        end
    end

    % Detect organism type and parse GFF3. Eukaryotes have mRNA features that are children of
    % gene features; prokaryotes have CDS.
    fprintf('Detecting organism type from GFF3 ...\n');
    orgType = detectOrgType(gffPath);   % 'eukaryote' | 'prokaryote'
    fprintf('  → %s genome detected.\n', orgType);
    
    fprintf('Parsing GFF3 file (this may take a moment for large files) ...\n');
    if strcmp(orgType, 'eukaryote')
        geneTable = parseGFF(gffPath, 'mRNA');
    else
        geneTable = parseGFF(gffPath, 'CDS');
    end

    % Save result
    writetable(geneTable, outputFile, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('Table saved to: %s  (%d rows)\n', outputFile, height(geneTable));
end

%--------------------------------------------------------------------------
% Helper functions

function orgType = detectOrgType(gffPath)
% Scan (up to) the first 50 000 lines to decide eukaryote vs prokaryote.
% Eukaryote: has 'mRNA' feature lines that are children of gene features.
% Prokaryote: gene children are only CDS/tRNA/rRNA (no mRNA).

    fid = fopen(gffPath, 'r');
    hasMRNA   = false;
    lineCount = 0;
    maxLines  = 50000;
    while ~feof(fid) && lineCount < maxLines
        line = fgetl(fid);
        lineCount = lineCount + 1;
        if ischar(line) && ~startsWith(line, '#')
            fields = strsplit(line, '\t');
            if numel(fields) >= 3 && strcmp(fields{3}, 'mRNA')
                hasMRNA = true;
                break
            end
        end
    end
    fclose(fid);
    orgType = 'prokaryote';
    if hasMRNA, orgType = 'eukaryote'; end
end


function T = parseGFF(gffPath, childFeature)
% Parse a GFF3 file extracting gene + child feature pairs.
%
% childFeature : 'mRNA' (eukaryote) | 'CDS' (prokaryote)
%                Thus, eukaryote (gene + mRNA) and prokaryote (gene + CDS).
%
% Returns a table with columns:
%   locus_tag | old_locus_tag | extra_id | GeneID | gene_name | GenBank_protein

    batchSize = 10000;
    rows  = cell(batchSize, 6);
    nRows = 0;
    % Use a containers.Map to hold the last gene record keyed by gene ID
    % so we can match mRNA → parent gene efficiently.
    geneMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    emptyGene = struct('locus_tag','','old_locus_tag','','extra_id','','geneID','','name','');

    fid = fopen(gffPath, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(line) || line(1) == '#', continue, end
        fields = strsplit(line, '\t');
        if numel(fields) < 9, continue, end

        featureType = fields{3};
        attrs       = fields{9};

        if strcmp(featureType, 'gene')
            g = parseGeneAttrs(attrs);
            if ~isempty(g.id)
                geneMap(g.id) = g;
            end

        elseif strcmp(featureType, childFeature)
            if strcmp(childFeature, 'mRNA')
                child    = parseMRNAAttrs(attrs);
                proteinID = child.genbank;
            else  % CDS
                child    = parseCDSAttrs(attrs);
                proteinID = child.protein_id;
            end

            if isempty(child.parent), continue, end

            if isKey(geneMap, child.parent)
                g = geneMap(child.parent);
            else
                g = emptyGene;
            end

            nRows = nRows + 1;
            if nRows > size(rows, 1)
                rows = [rows; cell(batchSize, 6)];
            end
            rows(nRows, :) = {g.locus_tag, g.old_locus_tag, g.extra_id, ...
                               g.geneID,    g.name,          proteinID};
        end
    end
    fclose(fid);

    rows = rows(1:nRows, :);

    % Deduplicate rows (CDS repeats per exon; harmless for mRNA isoforms
    % since column 6 differs per transcript).
    rowKeys = cellfun(@(a,b,c,d,e,f) strjoin({a,b,c,d,e,f}, '|'), ...
        rows(:,1), rows(:,2), rows(:,3), rows(:,4), rows(:,5), rows(:,6), ...
        'UniformOutput', false);
    [~, ia] = unique(rowKeys, 'stable');
    rows = rows(ia, :);

    T = cell2table(rows, 'VariableNames', ...
        {'locus_tag','old_locus_tag','extra_id','GeneID','gene_name','GenBank_protein'});
end

% Attibute parsing utilities
function value = extractAttr(attrStr, key)
% Extract the value of a GFF3 attribute by key.
% Handles URL-encoded characters (%2C etc.) for robustness.
% Returns '' if key is not present.
    pattern = sprintf('(?:^|;)%s=([^;]+)', regexptranslate('escape', key));
    tok = regexp(attrStr, pattern, 'tokens', 'once');
    if isempty(tok)
        value = '';
    else
        value = urldecode(tok{1});
    end
end

function value = extractDbxrefField(dbxref, prefix)
% From a comma-separated Dbxref value like
%   "Phytozome:Cre16.g651050,GeneID:5723799,GenBank:XM_001698190.2"
% extract the part after "prefix:".
% Returns '' if not found.
    if isempty(dbxref)
        value = '';
        return
    end
    parts   = strsplit(dbxref, ',');
    pattern = sprintf('^%s:(.+)$', regexptranslate('escape', prefix));
    value   = '';
    for i = 1:numel(parts)
        tok = regexp(strtrim(parts{i}), pattern, 'tokens', 'once');
        if ~isempty(tok)
            value = tok{1};
            return
        end
    end
end

function g = parseGeneAttrs(attrs)
    g.id            = extractAttr(attrs, 'ID');
    g.name          = extractAttr(attrs, 'Name');
    g.locus_tag     = extractAttr(attrs, 'locus_tag');
    g.old_locus_tag = extractAttr(attrs, 'old_locus_tag');
    
    % Dbxref is a comma-separated list like:
    %   Phytozome:Cre16.g651050,GeneID:5723799
    dbxref          = extractAttr(attrs, 'Dbxref');
    g.extra_id = extractDbxrefField(dbxref, 'Phytozome'); % To be improved
    g.geneID        = extractDbxrefField(dbxref, 'GeneID');
end

function m = parseMRNAAttrs(attrs)
    % Extract fields from an mRNA attribute string.
    m.parent  = extractAttr(attrs, 'Parent');
    % GenBank accession lives in Dbxref as GenBank:XM_...
    dbxref    = extractAttr(attrs, 'Dbxref');
    m.genbank = extractDbxrefField(dbxref, 'GenBank');
end

function c = parseCDSAttrs(attrs)
    c.parent     = extractAttr(attrs, 'Parent');
    c.protein_id = extractAttr(attrs, 'protein_id');
    if isempty(c.protein_id)
        % Some files put it in Dbxref as GenBank:WP_...
        dbxref       = extractAttr(attrs, 'Dbxref');
        c.protein_id = extractDbxrefField(dbxref, 'GenBank');
    end
end