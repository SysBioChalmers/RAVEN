function GSS = getUniProtScores(organism, varargin)
% getUniProtScores  Fetch curated subcellular locations from UniProt -> gene scoring structure.
%
% Queries the UniProtKB REST API for an organism's "Subcellular location" annotations and parses
% them into the GSS consumed by predictLocalization (no manual export needed). UniProt's curated
% location is qualitative, so each annotated compartment gets score 1.0 (a multi-located protein
% lands in several). Requires an internet connection.
%
% Parameters
% ----------
% organism : double or char
%     a UniProt organism / taxonomy id (e.g. 559292 for S. cerevisiae S288C).
%
% Name-Value Arguments
% --------------------
% idField : char
%     which identifier becomes the gene id: 'gene_oln' (ordered locus / ORF, the default, which
%     matches yeast-GEM), 'accession', or 'gene_primary'.
% reviewed : logical
%     restrict to curated Swiss-Prot entries (default true).
% compartmentMap : containers.Map
%     UniProt location term -> your model's compartment id (default
%     defaultCompartmentMap). Here the map is not only a renaming: its keys are
%     also the vocabulary searched for in UniProt's "Subcellular location [CC]"
%     free text, matched whole-word and case insensitively, so they must be
%     UniProt's own noun forms ('endoplasmic reticulum', not 'ER') and a term
%     that is not a key is never looked for. Several keys may share one
%     compartment id, which merges them; since every hit scores 1.0, merging
%     here means "annotated to any of these counts as this compartment". One
%     key cannot be split across two compartment ids - see parseScores for why,
%     and for the recipe to duplicate a column on the returned GSS instead.
% extraQuery : char
%     an additional UniProt query clause, ANDed into the query (e.g. 'gene:CIT1').
%
% Returns
% -------
% GSS : struct
%     a gene scoring structure for predictLocalization.
%
% Examples
% --------
%     GSS = getUniProtScores(559292);                       % all reviewed S. cerevisiae proteins
%     GSS = getUniProtScores(559292, 'idField', 'accession');
%     map = defaultCompartmentMap; map('cell wall') = 'ce';  % extend the vocabulary
%     GSS = getUniProtScores(559292, 'compartmentMap', map);
%
% See also
% --------
% parseScores, predictLocalization, defaultCompartmentMap

p = parseRAVENargs(varargin, {'idField','gene_oln'; 'reviewed',true; ...
                              'compartmentMap',[]; 'extraQuery',[]});
idField = char(p.idField);
if ~any(strcmp(idField, {'gene_oln','accession','gene_primary'}))
    error('RAVEN:badInput', 'idField must be gene_oln, accession, or gene_primary');
end

query = ['organism_id:' char(string(organism))];
if p.reviewed; query = [query ' AND reviewed:true']; end
if ~isempty(p.extraQuery); query = [query ' AND (' char(p.extraQuery) ')']; end

% Put the chosen id field first so parseScores reads it as the gene id column.
fields = [idField ',cc_subcellular_location'];
outFile = [tempname '.tsv'];
cleanup = onCleanup(@() i_deleteIfExists(outFile)); %#ok<NASGU> — removes temp file on return/error
websave(outFile, 'https://rest.uniprot.org/uniprotkb/stream', ...
        'query', query, 'format', 'tsv', 'fields', fields, weboptions('Timeout', 120));

% The UniProt /stream endpoint gzip-compresses the body and MATLAB does not auto-decompress it,
% so detect the gzip magic bytes and inflate before parsing.
fid = fopen(outFile, 'r');
magic = fread(fid, 2, 'uint8')';
fclose(fid);
if isequal(magic, [31 139])
    gz = [outFile '.gz'];
    movefile(outFile, gz);
    gunzip(gz);                 % writes outFile (drops the .gz extension)
    delete(gz);
end

map = p.compartmentMap;
if isempty(map); map = defaultCompartmentMap(); end
GSS = parseScores(outFile, 'predictor', 'uniprot', 'compartmentMap', map);
end

function i_deleteIfExists(f)
if isfile(f); delete(f); end
end
