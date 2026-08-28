function writeYAMLmodel(model,varargin)
% writeYAMLmodel  Write a model to a yaml file matching cobrapy's structure.
%
% The format is cobrapy's native !!omap layout, extended with RAVEN-only
% top-level per-entry keys (inchis, deltaG, metFrom, rxnFrom, references,
% confidence_score, protein) and the GECKO ec-rxns / ec-enzymes sections.
% Reaction EC numbers are written inside the `annotation` block as
% `ec-code` (the cobrapy/geckopy convention), not as a top-level reaction
% key. Output is byte-stable with raven_python's io.yaml.write_yaml_model
% when called with the same model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% fileName : char
%     name that the file will have. A dialog window will open if no file
%     name is specified.
% sortIds : logical
%     if metabolites, reactions, genes and compartments should be sorted
%     alphabetically by their identifier, otherwise they are kept in their
%     original order (default false).
%
% Examples
% --------
%     writeYAMLmodel(model,fileName,sortIds);
p=parseRAVENargs(varargin, {'fileName',[]; 'sortIds',false});
fileName=p.fileName; sortIds=p.sortIds;
if isempty(fileName)
    [fileName, pathName] = uiputfile({'*.yml;*.yaml'}, 'Select file for model export',[model.id '.yml']);
    if fileName == 0
        error('You should provide a file location')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);
if ~endsWith(fileName,{'.yml','.yaml'})
    fileName = strcat(fileName,'.yml');
end

%Check that model is in RAVEN format:
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Check that the model structure has no problems
checkModelStruct(model);

%Sort identifiers alphabetically
if sortIds == true
    model = sortIdentifiers(model);
end

%Simplify Miriam fields:
if isfield(model,'metMiriams')
    [model.newMetMiriams,model.newMetMiriamNames]   = extractMiriam(model.metMiriams);
end
if isfield(model,'rxnMiriams')
    [model.newRxnMiriams,model.newRxnMiriamNames]   = extractMiriam(model.rxnMiriams);
end
if isfield(model,'geneMiriams')
    [model.newGeneMiriams,model.newGeneMiriamNames] = extractMiriam(model.geneMiriams);
end
if isfield(model,'compMiriams')
    [model.newCompMiriams,model.newCompMiriamNames] = extractMiriam(model.compMiriams);
end

%Open file:
fid = fopen(fileName,'wt');
if fid == -1
    error(['Cannot write to ' fileName ', does the directory exist?'])
end
% cobrapy emits a bare `!!omap` root with no document-start marker;
% match that for byte-stable round-tripping.
fprintf(fid,'!!omap\n');

%Insert file header (metadata)
writeMetadata(model,fid);

%Metabolites:
% Field order matches cobrapy + raven_python.io.yaml:
%   id, name, compartment, charge, formula, annotation, then RAVEN-only
%   extras (inchis, deltaG, metFrom, notes).
% SMILES goes inside the annotation block (cobrapy convention), not at
% metabolite top level — the reader still accepts top-level `smiles:`
% for backward compatibility with older yeast-GEM files.
fprintf(fid,'- metabolites:\n');
for i = 1:length(model.mets)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'mets',        'txt', i, '    - id',          4)
    writeField(model, fid, 'metNames',    'txt', i, '    - name',        4)
    writeField(model, fid, 'metComps',    'txt', i, '    - compartment', 4)
    writeField(model, fid, 'metCharges',  'flt', i, '    - charge',      4)
    writeField(model, fid, 'metFormulas', 'txt', i, '    - formula',     4)
    writeAnnotation(model, fid, 'met',    i)
    writeField(model, fid, 'inchis',      'txt', i, '    - inchis',      4)
    writeField(model, fid, 'metDeltaG',   'flt', i, '    - deltaG',      4)
    writeField(model, fid, 'metFrom',     'txt', i, '    - metFrom',     4)
    writeField(model, fid, 'metNotes',    'txt', i, '    - notes',       4)
end

%Reactions:
% Field order matches cobrapy + raven_python.io.yaml:
%   id, name, metabolites, lower_bound, upper_bound, gene_reaction_rule,
%   objective_coefficient, subsystem, annotation (which carries the EC
%   numbers under `ec-code`, the cobrapy/geckopy convention), then
%   RAVEN-only extras (references, rxnFrom, deltaG, confidence_score,
%   notes). The notes key is the canonical `notes` (no longer
%   `rxnNotes`); the reader still accepts the legacy key.
fprintf(fid,'- reactions:\n');
for i = 1:length(model.rxns)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'rxns',                 'txt', i, '    - id',                    4)
    writeField(model, fid, 'rxnNames',             'txt', i, '    - name',                  4)
    writeField(model, fid, 'S',                    'txt', i, '    - metabolites',           4)
    writeField(model, fid, 'lb',                   'flt', i, '    - lower_bound',           4)
    writeField(model, fid, 'ub',                   'flt', i, '    - upper_bound',           4)
    writeField(model, fid, 'grRules',              'txtReq', i, '    - gene_reaction_rule', 4)
    if model.c(i)~=0
        writeField(model, fid, 'c',                'flt', i, '    - objective_coefficient', 4)
    end
    writeField(model, fid, 'subSystems',           'txt', i, '    - subsystem',             4)
    writeAnnotation(model, fid, 'rxn',             i)
    writeField(model, fid, 'rxnReferences',        'txt', i, '    - references',            4)
    writeField(model, fid, 'rxnFrom',              'txt', i, '    - rxnFrom',               4)
    writeField(model, fid, 'rxnDeltaG',            'flt', i, '    - deltaG',                4)
    writeField(model, fid, 'rxnConfidenceScores',  'flt', i, '    - confidence_score',      4)
    writeField(model, fid, 'rxnNotes',             'txt', i, '    - notes',                 4)
end

%Genes:
% genes is one of cobra's required top-level model keys (model_to_dict
% always emits it, empty or not) — write the flow-style empty list for
% a gene-less model, matching raven_toolbox.io.write_yaml_model exactly,
% rather than a bare `- genes:` with nothing following, which parses as
% `genes: null` and crashes readers that assume a present key means a
% list.
if isfield(model,'genes') && ~isempty(model.genes)
    fprintf(fid,'- genes:\n');
    for i = 1:length(model.genes)
        fprintf(fid,'  - !!omap\n');
        writeField(model, fid, 'genes',          'txt', i, '    - id',      4)
        writeField(model, fid, 'geneShortNames', 'txt', i, '    - name',    4)
        writeAnnotationSimple(model, fid, 'geneMiriams', 'newGeneMiriams', 'newGeneMiriamNames', i)
        writeField(model, fid, 'proteins',       'txt', i, '    - protein', 4)
    end
else
    fprintf(fid,'- genes: []\n');
end

%Compartments:
fprintf(fid,'- compartments: !!omap\n');
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames',   'txt', i, ['  - ' model.comps{i}], 2)
    writeAnnotationSimple(model, fid, 'compMiriams', 'newCompMiriams', 'newCompMiriamNames', i, '  ', 6)
end


%EC-model:
if isfield(model,'ec')
    % gecko_light flag at the top level (matches
    % raven_python.io.yaml — keeps the metaData block a pure provenance
    % container). The reader accepts both this key and the legacy
    % geckoLight key inside metaData.
    if model.ec.geckoLight; geckoLightStr = 'true'; else; geckoLightStr = 'false'; end
    fprintf(fid,'- gecko_light: %s\n', geckoLightStr);
    fprintf(fid,'- ec-rxns:\n');
    for i = 1:length(model.ec.rxns)
        fprintf(fid,'  - !!omap\n');
        writeField(model.ec, fid, 'rxns',      'txt', i, '    - id',      4)
        writeField(model.ec, fid, 'kcat',      'flt', i, '    - kcat',    4)
        writeField(model.ec, fid, 'source',    'txt', i, '    - source',  4)
        writeField(model.ec, fid, 'notes',     'txt', i, '    - notes',   4)
        writeField(model.ec, fid, 'eccodes',   'txt', i, '    - eccodes', 4)
        writeField(model.ec, fid, 'rxnEnzMat', 'txt', i, '    - enzymes', 4)
    end

    fprintf(fid,'- ec-enzymes:\n');
    for i = 1:length(model.ec.genes)
        fprintf(fid,'  - !!omap\n');
        writeField(model.ec, fid, 'genes',    'txt', i, '    - genes',    4)
        writeField(model.ec, fid, 'enzymes',  'txt', i, '    - enzymes',  4)
        writeField(model.ec, fid, 'mw',       'flt', i, '    - mw',       4)
        writeField(model.ec, fid, 'sequence', 'txt', i, '    - sequence', 4)
        writeField(model.ec, fid, 'concs',    'flt', i, '    - concs',    4)
    end
end

%Close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name,keyIndent)
%Writes a new line in the yaml file if the field exists and the field is
%not empty at the correspoinding position. It's recursive for the S
%field (reaction stoichiometry).

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end

    field = model.(fieldName);

    if strcmp(fieldName,'S')
        %S: create header & write each metabolite in a new line. Reactions
        %with no metabolites emit `metabolites: !!omap []` (the flow-style
        %empty omap cobrapy uses) so the file remains a valid YAML 1.2
        %document.
        if sum(field(:,pos) ~= 0) > 0
            fprintf(fid,'%s: !!omap\n',name);
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            %Sort metabolites:
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'flt', i, ['      - ' model.mets{i}], 6)
            end
        else
            fprintf(fid,'%s: !!omap []\n',name);
        end

    elseif strcmp(fieldName,'rxnEnzMat')
        %S: create header & write each enzyme in a new line
        fprintf(fid,'%s: !!omap\n',name);
        if sum(field(pos,:) ~= 0) > 0
            model.enzymes = model.enzymes(field(pos,:) ~= 0);
            model.coeffs  = field(pos,field(pos,:) ~= 0);
            %Sort metabolites:
            [model.enzymes,order] = sort(model.enzymes);
            model.coeffs          = model.coeffs(order);
            for i = 1:length(model.enzymes)
                writeField(model, fid, 'coeffs',  'flt', i, ['    - ' model.enzymes{i}], 4)
            end
        end

    elseif strcmp(fieldName,'subSystems')
        %subsystem: omitted entirely when the reaction carries none
        %(matching every other optional field), otherwise a block list —
        %always a list, even for a single subsystem, since a reaction can
        %have more than one and the reader already expects the list form.
        if iscell(field)
            if pos > numel(field); return; end
            list = field{pos};
        else
            if pos > size(field, 1); return; end
            list = field(pos, :);
        end
        if ischar(list); list = {list}; end
        list = list(~cellfun('isempty', list));
        if isempty(list)
            return
        end
        fprintf(fid,'%s:\n',name);
        indentStr = regexprep(name,'(^\s*).*','$1');
        for i = 1:length(list)
            emitListItem(fid, [indentStr '  -'], list{i}, keyIndent+2)
        end

    elseif strcmp(fieldName,'eccodes')
        %eccodes (GECKO ec-rxns section only): ;-joined string -> a
        %block list if more than one code, else a bare scalar.
        if isempty(field{pos})
            return
        end
        list = strip(strsplit(strrep(field{pos},' ',''),';'));
        list = list(~cellfun('isempty',list));
        if isempty(list)
            return
        elseif isscalar(list)
            emitScalarLine(fid, name, list{1}, keyIndent)
        else
            fprintf(fid,'%s:\n',name);
            indentStr = regexprep(name,'(^\s*).*','$1');
            for i = 1:length(list)
                emitListItem(fid, [indentStr '  -'], list{i}, keyIndent+2)
            end
        end

    elseif sum(pos) > 0
        %All other fields:
        if strcmp(type,'txt')
            value = field{pos};
            if ~isempty(value)
                emitScalarLine(fid, name, value, keyIndent)
            end
        elseif strcmp(type,'txtReq')
            %A structural field that is always emitted, even when empty
            %(matching cobra's "required reaction attribute" convention
            %for gene_reaction_rule) — unlike every optional field above,
            %which is simply omitted when empty.
            emitScalarLine(fid, name, field{pos}, keyIndent)
        elseif strcmp(type,'flt')
            if ~isnan(field(pos))
                emitScalarLine(fid, name, formatFloat(full(field(pos))), keyIndent, true)
            end
        end
    end
end
end

function writeMetadata(model, fid)
% Writes the metaData block. The `date` field is preserved when the model
% carries one (model.date), so round-trips don't churn on every write; if
% absent it's filled with the current date.

fprintf(fid, '- metaData: !!omap\n');
emitMetaField(fid, 'id',           valueOrDefault(model,'id','blankID'));
emitMetaField(fid, 'name',         valueOrDefault(model,'name','blankName'));
if isfield(model,'version') && ~isempty(model.version)
    emitMetaField(fid, 'version', model.version);
end
if isfield(model,'date') && ~isempty(model.date)
    dateValue = model.date;
else
    dateValue = datestr(now, 29); %#ok<DATST>  % 29 = yyyy-mm-dd
end
emitMetaField(fid, 'date', dateValue);
if isfield(model,'annotation')
    annoFields = {'defaultLB','defaultUB','givenName','familyName', ...
                  'authors','email','organization','taxonomy','note','sourceUrl'};
    for k = 1:numel(annoFields)
        f = annoFields{k};
        if isfield(model.annotation, f) && ~isempty(model.annotation.(f))
            emitMetaField(fid, f, model.annotation.(f));
        end
    end
end
% gecko_light is emitted at the top level (see geckoLight emission near
% the GECKO ec-* sections) to match raven_python.io.yaml; keeping it
% out of metaData lets cobrapy/ruamel keep the section a pure
% provenance block.
end

function v = valueOrDefault(model, field, defaultVal)
if isfield(model, field) && ~isempty(model.(field))
    v = model.(field);
else
    v = defaultVal;
end
end

function emitMetaField(fid, key, value)
% Emit one `  - key: value` line inside the metaData omap block.
if islogical(value)
    if value; value = 'true'; else; value = 'false'; end
    emitScalarLine(fid, ['  - ' key], value, 2, true);
elseif isnumeric(value)
    % defaultLB/defaultUB are genuinely floats (cobra's Reaction bounds);
    % every other metaData annotation value that reaches here numeric
    % (e.g. a stray numeric taxonomy id) is formatted the same way, since
    % none of RAVEN's own metaData fields are conceptually integers.
    emitScalarLine(fid, ['  - ' key], formatFloat(double(value)), 2, true);
else
    if ~ischar(value) && ~isstring(value)
        value = char(value);
    end
    emitScalarLine(fid, ['  - ' key], char(value), 2);
end
end

function s = formatFloat(x)
% Format every numeric field the same way: a whole number still gets an
% explicit ".0" (e.g. "1000.0", not "1000"; "2.0" for a charge, not "2").
if isinf(x)
    if x > 0
        s = '.inf';
    else
        s = '-.inf';
    end
    return;
end
if x == floor(x) && abs(x) < 1e16
    s = sprintf('%.1f', x);   % e.g. 1000 -> "1000.0"
else
    s = sprintf('%.15g', x);
end
end

function writeAnnotationSimple(model, fid, miriamsField, miriamValsField, miriamNamesField, pos, indent, subIndent)
% Emit a bare `annotation: !!omap` block from MIRIAM cross-references only
% (used for genes and compartments, which carry no extra non-MIRIAM
% annotation key; cobrapy has no equivalent for either — compartments in
% particular are a flat {id: name} map on the cobra side, with no
% per-compartment annotation at all — so there is no reference format to
% match here, only readYAMLmodel's own expectation that a compartment's
% annotation sub-entries are indented >=6 spaces). Keys are emitted in
% alphabetical order, matching cobra's sorted-dict annotation elsewhere.
if nargin < 7
    indent = '    ';
end
if nargin < 8
    subIndent = length(indent)+2;
end
if ~isfield(model, miriamsField) || isempty(model.(miriamsField){pos})
    return
end
names  = model.(miriamNamesField);
values = model.(miriamValsField);
entries = collectMiriamEntries(names, values, pos);
if isempty(entries)
    return
end
[~,order] = sort({entries.name});
entries = entries(order);
fprintf(fid, '%s- annotation: !!omap\n', indent);
emitAnnotationEntries(fid, entries, subIndent);
end

function entries = collectMiriamEntries(names, values, pos)
entries = struct('name', {}, 'list', {}, 'forceList', {});
for j = 1:size(values, 2)
    v = values{pos, j};
    if isempty(v); continue; end
    list = strip(strsplit(strrep(v, ' ', ''), ';'));
    entries(end+1) = struct('name', names{j}, 'list', {list}, 'forceList', false); %#ok<AGROW>
end
end

function writeAnnotation(model, fid, kind, pos)
% Emit the per-entry `annotation` block, fusing MIRIAM cross-references
% with the non-MIRIAM cobrapy-style annotation keys: `smiles` for
% metabolites and `ec-code` (EC numbers) for reactions. cobrapy and
% geckopy read both from inside `annotation` (not as top-level entry
% keys), so this helper keeps the YAML aligned. `ec-code` is emitted as a
% list — matching cobrapy / raven-python and geckopy, which read it as
% list[str] — even when there is a single code. All keys (MIRIAMs and the
% extra key together) are emitted in alphabetical order, matching cobra's
% sorted-dict annotation.
switch kind
    case 'met'
        miriamsField     = 'metMiriams';
        miriamNamesField = 'newMetMiriamNames';
        miriamValsField  = 'newMetMiriams';
        extraName        = 'smiles';
        extraField       = 'metSmiles';
    case 'rxn'
        miriamsField     = 'rxnMiriams';
        miriamNamesField = 'newRxnMiriamNames';
        miriamValsField  = 'newRxnMiriams';
        extraName        = 'ec-code';
        extraField       = 'eccodes';
    otherwise
        error('writeAnnotation:unsupportedKind', 'Unsupported kind: %s', kind);
end

hasMiriams = isfield(model, miriamsField) && ~isempty(model.(miriamsField){pos});
hasExtra   = isfield(model, extraField) && pos <= numel(model.(extraField)) ...
             && ~isempty(model.(extraField){pos});

if ~hasMiriams && ~hasExtra
    return;
end

entries = struct('name', {}, 'list', {}, 'forceList', {});
if hasMiriams
    entries = collectMiriamEntries(model.(miriamNamesField), model.(miriamValsField), pos);
end
if hasExtra
    % EC numbers / SMILES are always emitted as a list (matching
    % cobrapy/raven-python/geckopy, which read `ec-code` as list[str]),
    % even when there is a single value — unlike a MIRIAM entry, which
    % collapses to a bare scalar when there is only one.
    codes = strip(strsplit(strrep(model.(extraField){pos}, ' ', ''), ';'));
    codes = codes(~cellfun('isempty', codes));
    if ~isempty(codes)
        entries(end+1) = struct('name', extraName, 'list', {codes}, 'forceList', true); %#ok<AGROW>
    end
end
if isempty(entries)
    return
end
[~,order] = sort({entries.name});
entries = entries(order);

fprintf(fid, '    - annotation: !!omap\n');
emitAnnotationEntries(fid, entries, 6);
end

function emitAnnotationEntries(fid, entries, keyIndent)
% Shared tail for writeAnnotation / writeAnnotationSimple: entries are
% already sorted by name; emit each as a bare scalar (single value) or a
% block list (multiple values), at the given key indent.
pad = repmat(' ', 1, keyIndent);
for e = entries
    if numel(e.list) == 1 && ~e.forceList
        emitScalarLine(fid, [pad '- ' e.name], e.list{1}, keyIndent)
    else
        fprintf(fid, '%s- %s:\n', pad, e.name);
        for k = 1:numel(e.list)
            emitListItem(fid, [pad '  -'], e.list{k}, keyIndent+2)
        end
    end
end
end

function emitScalarLine(fid, prefix, value, keyIndent, bareNumeric)
% Emit "PREFIX: VALUE\n" (double-quoted when needed), where PREFIX is the
% "    - key" part with no trailing colon and value is the raw,
% unescaped scalar content.
%
% bareNumeric (default false): value is already a pre-formatted number
% (from formatFloat) that must never be quoted — numbers are never
% ambiguous once formatted.
if nargin < 5
    bareNumeric = false;
end
fprintf(fid, '%s:', prefix);
emitScalarCore(fid, length(prefix)+1, value, keyIndent, bareNumeric)
end

function emitListItem(fid, dashPrefix, value, keyIndent)
% Emit "DASHPREFIX VALUE\n" for a plain sequence entry (no key, just the
% "- " marker) — e.g. a subsystem, EC code, SMILES or MIRIAM list item.
% Unlike emitScalarLine, no colon is written.
fprintf(fid, '%s', dashPrefix);
emitScalarCore(fid, length(dashPrefix), value, keyIndent, false)
end

function emitScalarCore(fid, ~, value, ~, bareNumeric)
% Shared tail for emitScalarLine / emitListItem: writes a single space
% (as ruamel always does after ':' or '-'), then the value, double-quoted
% when needed (matching Prettier's YAML default, rather than ruamel's own
% single-quote default). No folding — RAVEN/raven-toolbox's shared format
% never wraps a scalar across lines, however long, so the two unused
% positional arguments (column and key indent) only exist to keep every
% call site's signature stable.
if bareNumeric
    fprintf(fid, ' %s\n', value);
    return
end
if needsQuote(value)
    escaped = strrep(value, '\', '\\');
    escaped = strrep(escaped, '"', '\"');
    fprintf(fid, ' "%s"\n', escaped);
else
    fprintf(fid, ' %s\n', value);
end
end

function tf = needsQuote(s)
% Whether ruamel would quote this scalar: either because a plain (bare)
% reading would resolve to a non-string YAML type (bool/int/float/null/
% timestamp — cobra's dumper never emits explicit type tags), or because
% the text itself is not valid as a plain block scalar (leading
% indicator character, a mid-string ": "/" #", leading/trailing
% whitespace, or embedded control characters). Ports the relevant
% subset of ruamel.yaml's Resolver.implicit_resolvers and
% Emitter.analyze_scalar (YAML 1.2 variants; cobra's round-trip dumper
% does not pin an older version).
if isempty(s)
    tf = true;
    return
end

resolverPattern = [ ...
    '^(true|True|TRUE|false|False|FALSE)$|' ...
    '^[-+]?[0-9][0-9_]*\.[0-9_]*([eE][-+]?[0-9]+)?$|' ...
    '^[-+]?[0-9][0-9_]*[eE][-+]?[0-9]+$|' ...
    '^[-+]?\.[0-9_]+([eE][-+][0-9]+)?$|' ...
    '^[-+]?\.(inf|Inf|INF)$|' ...
    '^\.(nan|NaN|NAN)$|' ...
    '^[-+]?0b[01_]+$|' ...
    '^[-+]?0o?[0-7_]+$|' ...
    '^[-+]?[0-9_]+$|' ...
    '^[-+]?0x[0-9a-fA-F_]+$|' ...
    '^(~|[Nn]ull|NULL)$|' ...
    '^[0-9]{4}-[0-9]{2}-[0-9]{2}$' ...
];
if ~isempty(regexp(s, resolverPattern, 'once'))
    tf = true;
    return
end

first = s(1);
if any(first == '#,[]{}&*!|>''"%@`')
    tf = true; return
end
followedByWs = (numel(s) == 1) || isWs(s(2));
if (first == '?' || first == ':') && followedByWs
    tf = true; return
end
if first == '-' && followedByWs
    tf = true; return
end
if isWs(first) || isWs(s(end))
    tf = true; return
end
if any(s == char(10)) || any(s == char(13))
    tf = true; return
end
for i = 2:numel(s)
    ch = s(i);
    if ch == ':' && ((i == numel(s)) || isWs(s(i+1)))
        tf = true; return
    end
    if ch == '#' && isWs(s(i-1))
        tf = true; return
    end
end
tf = false;
end

function tf = isWs(ch)
tf = ch == ' ' || ch == char(9) || ch == char(10) || ch == char(13);
end
