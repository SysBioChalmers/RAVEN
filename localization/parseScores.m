function GSS = parseScores(inputFile, varargin)
% parseScores  Parse predictor / database output into a gene scoring structure (GSS).
%
% Reads the output of a subcellular-localization predictor or evidence database and
% builds the gene scoring structure consumed by predictLocalization. Scores are
% normalized so the best score for each gene is 1.0.
%
% Parameters
% ----------
% inputFile : char
%     the file with the predictor / database output.
%
% Name-Value Arguments
% --------------------
% predictor : char
%     'deeploc' (DeepLoc 2 per-protein CSV, default), 'cello' (CELLO),
%     'mulocdeep' (MULocDeep wide table: gene id + one probability column per
%     compartment), 'compartments' (a COMPARTMENTS jensenlab.org channel TSV), or
%     'uniprot' (a UniProtKB TSV export with a "Subcellular location [CC]" column).
% compartmentMap : containers.Map
%     maps predictor/database compartment labels to your model's compartment ids and
%     merges synonyms (max). For 'uniprot' it is also the vocabulary of location terms
%     to search for and defaults to defaultCompartmentMap; for the others it is optional
%     (labels not in the map are dropped when one is given).
% idColumn : char or double
%     column (name or 1-based index) holding the gene id for 'mulocdeep' / 'uniprot'.
%     Default: the first column. For yeast-GEM, use the UniProt ordered-locus column.
% minConfidence : double
%     for 'compartments', drop annotations with confidence below this (default 0).
%
% Returns
% -------
% GSS : struct
%     a gene scoring structure (genes, compartments, scores) for predictLocalization.
%
% Examples
% --------
%     GSS = parseScores(file, 'predictor', 'deeploc');
%     GSS = parseScores(file, 'predictor', 'uniprot', 'idColumn', 'Gene Names (ordered locus)');
%
% See also
% --------
% predictLocalization, getUniProtScores, defaultCompartmentMap

p = parseRAVENargs(varargin, {'predictor',[]; 'compartmentMap',[]; ...
                              'idColumn',[]; 'minConfidence',0});
predictor = p.predictor;
if isempty(predictor); predictor = 'deeploc'; else; predictor = char(predictor); end
compartmentMap = p.compartmentMap;

if ~isfile(inputFile)
    error('RAVEN:badInput', 'Could not open file %s', char(string(inputFile)));
end

switch lower(predictor)
    case 'cello'
        GSS = i_parseCello(inputFile);
    case 'deeploc'
        GSS = i_parseDeepLoc(inputFile);
    case 'mulocdeep'
        GSS = i_parseWide(inputFile, p.idColumn);
    case 'compartments'
        GSS = i_parseCompartments(inputFile, p.minConfidence);
    case 'uniprot'
        if isempty(compartmentMap); compartmentMap = defaultCompartmentMap(); end
        GSS = i_parseUniProt(inputFile, p.idColumn, compartmentMap);
        compartmentMap = [];   % already applied during parsing
    case 'wolf'
        error('RAVEN:badInput', ['WoLF PSORT support has been removed. Use a modern ' ...
            'predictor (''deeploc'', ''mulocdeep''), the ''compartments'' database, or ' ...
            '''uniprot'' (see getUniProtScores).']);
    otherwise
        error('RAVEN:badInput', 'Unknown predictor "%s"', predictor);
end

% Remove duplicate genes
[~, J, K] = unique(GSS.genes);
if numel(J) ~= numel(K)
    warning('RAVEN:warning', 'There are duplicate genes in the input file');
    GSS.genes = GSS.genes(J);
    GSS.scores = GSS.scores(J,:);
end

% Optional renaming of compartments to model ids (already done for 'uniprot')
if ~isempty(compartmentMap)
    GSS = i_applyCompartmentMap(GSS, compartmentMap);
end

if isempty(GSS.genes) || isempty(GSS.compartments)
    warning('RAVEN:warning', ['parseScores produced an empty gene scoring structure ' ...
        '(no genes or no compartments matched). Check the input file and, if used, that ' ...
        'compartmentMap labels match the predictor/database compartment names.']);
end

% Normalize so the best score per gene is 1.0 (genes with no positive evidence untouched)
m = max(GSS.scores, [], 2);
m(m <= 0) = 1;
GSS.scores = bsxfun(@times, GSS.scores, 1./m);
end

% ------------------------------------------------------------------------- parsers

function GSS = i_parseCello(inputFile)
fid = fopen(inputFile, 'r');
tline = fgetl(fid);
tline = regexprep(tline, '^.+#Combined:\t', '');
tline = regexprep(tline, '\t#Most-likely-Location.+', '');
GSS.compartments = transpose(regexp(tline, '\t', 'split'));
GSS.genes = {}; GSS.scores = [];
row = 0;
while true
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    tline = regexprep(tline, '^.+:\t', '');
    tline = regexprep(tline, ' .+', '');
    tline = regexp(tline, '\t', 'split');
    if numel(tline) < numel(GSS.compartments)
        error('RAVEN:badInput', 'Malformed CELLO row %d (fewer fields than compartments)', row+1);
    end
    row = row + 1;
    GSS.scores(row,:) = str2double(tline(1:numel(GSS.compartments)));
    GSS.genes{row,1} = tline{1,end};
end
fclose(fid);
end

function GSS = i_parseDeepLoc(inputFile)
fid = fopen(inputFile, 'r');
tline = fgetl(fid);
comps = regexp(tline, ',', 'split');
GSS.compartments = transpose(comps(4:end));   % cols 1-3 are Protein_ID/Localizations/Signals
GSS.genes = {}; GSS.scores = [];
row = 0;
while true
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    parts = regexp(tline, ',', 'split');
    row = row + 1;
    GSS.scores(row,:) = str2double(parts(4:end));
    GSS.genes{row,1} = parts{1,1};
end
fclose(fid);
end

function GSS = i_parseWide(inputFile, idColumn)
% MULocDeep / generic wide table: gene id column + one numeric column per compartment.
T = readtable(inputFile, 'FileType', 'text', 'VariableNamingRule', 'preserve');
idx = i_idColumnIndex(T, idColumn);
GSS.genes = cellstr(string(T{:,idx}));
rest = setdiff(1:width(T), idx, 'stable');
GSS.compartments = transpose(T.Properties.VariableNames(rest));
GSS.scores = i_numericMatrix(T, rest);
end

function GSS = i_parseCompartments(inputFile, minConfidence)
% COMPARTMENTS channel TSV (no header): col1 gene, col4 GO term name, last col confidence.
T = readtable(inputFile, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
if width(T) < 5
    error('RAVEN:badInput', 'Expected a COMPARTMENTS full TSV (>=5 columns)');
end
genes = string(T{:,1});
comps = string(T{:,4});
score = T{:,end};
if ~isnumeric(score); score = str2double(string(score)); end
keep = ~isnan(score) & score >= minConfidence;
genes = genes(keep); comps = comps(keep); score = score(keep);

ug = unique(genes, 'stable');
uc = unique(comps, 'stable');
M = zeros(numel(ug), numel(uc));
[~, gi] = ismember(genes, ug);
[~, ci] = ismember(comps, uc);
for k = 1:numel(genes)
    if M(gi(k), ci(k)) < score(k); M(gi(k), ci(k)) = score(k); end
end
GSS.genes = cellstr(ug);
GSS.compartments = cellstr(uc);
GSS.scores = M;
end

function GSS = i_parseUniProt(inputFile, idColumn, map)
% UniProtKB TSV export: scan the "Subcellular location [CC]" text for known terms.
T = readtable(inputFile, 'FileType', 'text', 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
vn = string(T.Properties.VariableNames);
locIdx = find(contains(lower(vn), 'subcellular'), 1);
if isempty(locIdx)
    error('RAVEN:badInput', 'No "Subcellular location" column in the UniProt file');
end
idIdx = i_idColumnIndex(T, idColumn);
ids = string(T{:,idIdx});
loc = string(T{:,locIdx});

labels = lower(string(keys(map)));
codesForLabels = string(values(map));
codeList = unique(codesForLabels, 'stable');
% Whole-word patterns so a short label (e.g. 'cytoplasm') does not match inside a longer,
% distinct location term (e.g. 'Cytoplasmic vesicle'). Labels must be the UniProt noun forms.
patterns = "\<" + arrayfun(@(s) string(regexptranslate('escape', char(s))), labels) + "\>";

M = zeros(numel(ids), numel(codeList));
hasHit = false(numel(ids), 1);
for i = 1:numel(ids)
    txt = lower(loc(i));
    txt = regexprep(txt, '\{[^}]*\}', ' ');         % drop {ECO:...} evidence
    % Split into per-isoform blocks and drop the Note=... free text WITHIN each block, so a
    % later isoform's locations survive and Note prose is never read as a location.
    blocks = regexp(char(txt), '\[isoform[^\]]*\]:', 'split');
    for b = 1:numel(blocks)
        seg = regexprep(blocks{b}, 'note=.*', ' ');
        for j = 1:numel(labels)
            if ~isempty(regexp(seg, patterns(j), 'once'))
                c = find(codeList == codesForLabels(j), 1);
                M(i, c) = 1; hasHit(i) = true;
            end
        end
    end
end
keep = hasHit & strlength(strtrim(ids)) > 0;
GSS.genes = cellstr(ids(keep));
GSS.compartments = cellstr(codeList);
GSS.scores = M(keep, :);
end

% ------------------------------------------------------------------------- helpers

function idx = i_idColumnIndex(T, idColumn)
if isempty(idColumn)
    idx = 1;
elseif isnumeric(idColumn)
    idx = idColumn;
else
    idx = find(strcmp(T.Properties.VariableNames, char(idColumn)), 1);
    if isempty(idx)
        error('RAVEN:badInput', 'idColumn "%s" not found', char(idColumn));
    end
end
end

function M = i_numericMatrix(T, cols)
M = zeros(height(T), numel(cols));
for k = 1:numel(cols)
    col = T{:, cols(k)};
    if ~isnumeric(col); col = str2double(string(col)); end
    M(:, k) = col;
end
end

function GSS = i_applyCompartmentMap(GSS, map)
labels = lower(string(keys(map)));
codes = string(values(map));
old = lower(string(GSS.compartments));
mapped = strings(numel(old), 1);
keepCol = false(numel(old), 1);
for i = 1:numel(old)
    j = find(labels == old(i), 1);
    if ~isempty(j); mapped(i) = codes(j); keepCol(i) = true; end
end
usedCodes = unique(mapped(keepCol), 'stable');
M = zeros(numel(GSS.genes), numel(usedCodes));
for c = 1:numel(usedCodes)
    cols = keepCol & mapped == usedCodes(c);
    M(:, c) = max(GSS.scores(:, cols), [], 2);
end
GSS.compartments = cellstr(usedCodes);
GSS.scores = M;
end
