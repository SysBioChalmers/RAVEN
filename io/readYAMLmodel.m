function model=readYAMLmodel(fileName, verbose)
% readYAMLmodel
%   Reads a model from a YAML file. Two on-disk shapes are accepted:
%
%   1. The canonical geckopy schema (docs/yaml_format.md in geckopy):
%      top-level mapping with `id`, `name`, optional `metaData`,
%      optional `gecko_light`, `compartments` as a flat mapping,
%      `metabolites`/`reactions`/`genes` as lists of plain mappings,
%      reaction `metabolites` and ec-rxn `enzymes` as flat mappings,
%      and `annotation` blocks as namespace-keyed lists of strings.
%      This is what the current writeYAMLmodel emits.
%
%   2. The legacy RAVEN format (outer `---`/`!!omap` wrapper, `!!omap`
%      tags throughout, `id`/`name`/`geckoLight` nested under
%      `metaData`, scalar annotation values). Older YAML files written
%      by past RAVEN/GECKO releases use this shape.
%
%   Format detection is automatic; both shapes yield the same model
%   struct on the MATLAB side.
%
%   Input:
%   fileName    a model file in yaml file format. A dialog window will open
%               if no file name is specified.
%   verbose     set as true to monitor progress (optional, default false)
%
%   Output:
%   model       a model structure
%
% Usage: model = readYAMLmodel(fileName, verbose)
if nargin<1 || isempty(fileName)
    [fileName, pathName] = uigetfile({'*.yml;*.yaml'}, 'Please select the model file');
    if fileName == 0
        error('You should select a model file')
    else
        fileName = fullfile(pathName,fileName);
    end
end
if nargin < 2
    verbose = false;
end

if ~isfile(fileName)
    error('Yaml file %s cannot be found', string(fileName));
end

if verLessThan('matlab','9.9') %readlines introduced 2020b
    fid=fopen(fileName);
    line_raw=cell(1000000,1);
    i=1;
    while ~feof(fid)
        line_raw{i}=fgetl(fid);
        i=i+1;
    end
    line_raw(i:end)=[];
    line_raw=string(line_raw);
    fclose(fid);
else
    line_raw=readlines(fileName);
end

% Dispatch on detected format. Legacy files open with `---` and use
% `!!omap` tags; canonical files do not.
if isLegacyFormat(line_raw)
    model = parseLegacy(line_raw, verbose);
else
    model = parseCanonical(line_raw, verbose);
end
end

% --- Format detection ---------------------------------------------------

function tf = isLegacyFormat(line_raw)
tf = false;
for k = 1:min(numel(line_raw), 50)
    s = strtrim(char(line_raw{k}));
    if isempty(s)
        continue
    end
    if strcmp(s,'---') || strcmp(s,'!!omap') || contains(s,'!!omap')
        tf = true;
        return
    end
    % First non-empty, non-comment line that's neither `---` nor `!!omap`
    % means we're in canonical land.
    if s(1) ~= '#'
        return
    end
end
end

% --- Canonical-format parser -------------------------------------------

function model = parseCanonical(line_raw, verbose)
%Indent-based parser for the canonical geckopy YAML schema. Strategy: a
%small recursive-descent parser that walks the line stream, tracking the
%current indentation context. We materialise the YAML as a nested
%struct/cell tree, then translate it to a RAVEN model struct.

if verbose
    fprintf('Parsing YAML (canonical schema)...\n');
end

%Strip blank lines and full-line comments so the parser only sees real
%content. We keep the line index 1-based.
lines = string(line_raw);
keep  = ~arrayfun(@(s) isempty(strtrim(char(s))) || startsWith(strtrim(char(s)),'#'), lines);
lines = lines(keep);

[doc, ~] = parseMapping(cellstr(lines), 1, 0);

model = canonicalToModel(doc, verbose);
end

function [out, idx] = parseMapping(lines, idx, baseIndent)
%Parse a YAML mapping whose keys live at column == baseIndent. Stops as
%soon as we hit a line with less indentation than baseIndent or the EOF.
out = struct();
while idx <= numel(lines)
    raw = lines{idx};
    [indent, content] = splitIndent(raw);
    if indent < baseIndent
        return
    end
    if indent > baseIndent
        %Defensive: a mapping line at unexpected indent. Skip it.
        idx = idx + 1;
        continue
    end
    if startsWith(content,'-')
        %A list under a mapping key would have been consumed by the parent
        %call when it saw `key:` with no inline value. Reaching `-` here
        %means we've stepped outside the mapping; return.
        return
    end
    [key, value, hasInline] = splitKeyValue(content);
    if isempty(key)
        idx = idx + 1;
        continue
    end
    safeKey = sanitiseKey(key);
    if hasInline
        out.(safeKey) = parseScalar(value);
        idx = idx + 1;
    else
        %Look ahead at the value block. Cases:
        %  next line indent > baseIndent, starts with `-` => nested list
        %  next line indent > baseIndent, no `-` => nested mapping
        %  next line indent == baseIndent, starts with `-` => "compact"
        %      list (YAML lets a list share indent with its parent key,
        %      the standard cobra-style annotation idiom)
        %  anything else => empty value
        nextIdx = idx + 1;
        if nextIdx > numel(lines)
            out.(safeKey) = '';
            return
        end
        [nIndent, nContent] = splitIndent(lines{nextIdx});
        if nIndent == baseIndent && startsWith(nContent,'-')
            [child, idx] = parseList(lines, nextIdx, baseIndent);
            out.(safeKey) = child;
        elseif nIndent > baseIndent
            if startsWith(nContent,'-')
                [child, idx] = parseList(lines, nextIdx, nIndent);
            else
                [child, idx] = parseMapping(lines, nextIdx, nIndent);
            end
            out.(safeKey) = child;
        else
            out.(safeKey) = '';
            idx = nextIdx;
        end
    end
end
end

function [out, idx] = parseList(lines, idx, baseIndent)
%Parse a YAML list whose `-` markers live at column == baseIndent.
out = {};
while idx <= numel(lines)
    raw = lines{idx};
    [indent, content] = splitIndent(raw);
    if indent < baseIndent
        return
    end
    if indent > baseIndent
        idx = idx + 1;
        continue
    end
    if ~startsWith(content,'-')
        return
    end
    %Strip the `- ` marker, leaving either a scalar or the first
    %key/value of a mapping element.
    after = regexprep(content,'^-\s*','');
    if isempty(after)
        %Item is a bare `-` followed by a nested mapping on the next line.
        nextIdx = idx + 1;
        if nextIdx > numel(lines)
            out{end+1,1} = ''; %#ok<AGROW>
            return
        end
        [nIndent, nContent] = splitIndent(lines{nextIdx});
        if startsWith(nContent,'-')
            [child, idx] = parseList(lines, nextIdx, nIndent);
        else
            [child, idx] = parseMapping(lines, nextIdx, nIndent);
        end
        out{end+1,1} = child; %#ok<AGROW>
    else
        [key, value, hasInline] = splitKeyValue(after);
        if isempty(key)
            %Scalar list item, e.g. `- foo`.
            out{end+1,1} = parseScalar(after); %#ok<AGROW>
            idx = idx + 1;
        else
            %The element is a mapping; its first key sits inline with the
            %`-`. Capture this first key, then if there are more keys
            %indented further, parse them as the remainder of the mapping.
            element = struct();
            safeKey = sanitiseKey(key);
            if hasInline
                element.(safeKey) = parseScalar(value);
                idx = idx + 1;
                %Continuation: subsequent lines indented at (baseIndent + 2)
                %belong to the same element.
                while idx <= numel(lines)
                    raw2 = lines{idx};
                    [i2, c2] = splitIndent(raw2);
                    if i2 < baseIndent + 2
                        break
                    end
                    if i2 > baseIndent + 2
                        idx = idx + 1;
                        continue
                    end
                    if startsWith(c2,'-')
                        %Reached the next list element belonging to a
                        %parent list (e.g. annotation values). Stop here.
                        break
                    end
                    [k2, v2, has2] = splitKeyValue(c2);
                    if isempty(k2)
                        idx = idx + 1;
                        continue
                    end
                    sk2 = sanitiseKey(k2);
                    if has2
                        element.(sk2) = parseScalar(v2);
                        idx = idx + 1;
                    else
                        nxt = idx + 1;
                        if nxt > numel(lines)
                            element.(sk2) = '';
                            idx = nxt;
                            break
                        end
                        [ni, nc] = splitIndent(lines{nxt});
                        if ni <= i2
                            element.(sk2) = '';
                            idx = nxt;
                            continue
                        end
                        if startsWith(nc,'-')
                            [grand, idx] = parseList(lines, nxt, ni);
                        else
                            [grand, idx] = parseMapping(lines, nxt, ni);
                        end
                        element.(sk2) = grand;
                    end
                end
            else
                %Inline key with no value -> nested block on following lines.
                nextIdx = idx + 1;
                if nextIdx > numel(lines)
                    element.(safeKey) = '';
                    idx = nextIdx;
                else
                    [nIndent, nContent] = splitIndent(lines{nextIdx});
                    if startsWith(nContent,'-')
                        [grand, idx] = parseList(lines, nextIdx, nIndent);
                    else
                        [grand, idx] = parseMapping(lines, nextIdx, nIndent);
                    end
                    element.(safeKey) = grand;
                end
            end
            out{end+1,1} = element; %#ok<AGROW>
        end
    end
end
end

function [indent, content] = splitIndent(line)
line = char(line);
n = numel(line);
i = 1;
while i <= n && line(i) == ' '
    i = i + 1;
end
indent = i - 1;
content = line(i:end);
end

function [key, value, hasInline] = splitKeyValue(content)
%Detect `key: value` (hasInline=true), `key:` (hasInline=false), or no key
%(returns empty key). A `:` inside a quoted string is not a separator.
key = '';
value = '';
hasInline = false;
content = char(content);
if isempty(content)
    return
end
%Find the first colon that's not inside quotes.
inQuote = false;
colonPos = 0;
for k = 1:numel(content)
    c = content(k);
    if c == '"'
        inQuote = ~inQuote;
    elseif c == ':' && ~inQuote
        colonPos = k;
        break
    end
end
if colonPos == 0
    return
end
key = strtrim(content(1:colonPos-1));
rest = strtrim(content(colonPos+1:end));
if isempty(rest)
    hasInline = false;
else
    hasInline = true;
    value = rest;
end
end

function v = parseScalar(s)
s = strtrim(char(s));
if isempty(s)
    v = '';
    return
end
%Strip an inline `# ...` comment (only when `#` is preceded by whitespace
%and we're not inside a quoted string).
if s(1) ~= '"'
    commentStart = regexp(s,'\s#','once');
    if ~isempty(commentStart)
        s = strtrim(s(1:commentStart));
    end
end
%Strip surrounding double quotes.
if numel(s) >= 2 && s(1) == '"' && s(end) == '"'
    v = s(2:end-1);
    return
end
%Normalise YAML's `.nan`/`.NaN`/`nan` to the canonical 'NaN' literal so
%downstream str2double calls produce NaN.
if strcmpi(s,'.nan') || strcmpi(s,'nan')
    v = 'NaN';
    return
end
v = s;
end

function s = sanitiseKey(key)
%Encode YAML keys as legal MATLAB field names. We escape just `.` and
%`-` (the two character classes that show up in cobra-style namespace
%keys like `bigg.metabolite` and `ec-code`) and prefix a leading-digit
%guard so BiGG-style ids like `2agpe120` round-trip cleanly. Literal
%`_dot_`/`_dash_`/`kbeg_` substrings in original keys are not protected
%against; they're vanishingly rare in practice.
s = strrep(key,'.','_dot_');
s = strrep(s,'-','_dash_');
s = regexprep(s,'[^A-Za-z0-9_]','_');
if isempty(s) || ~isstrprop(s(1),'alpha')
    s = ['kbeg_' s];
end
end

% --- Canonical YAML tree -> RAVEN model struct -------------------------

function model = canonicalToModel(doc, verbose)
%Translate the intermediate struct/cell tree produced by parseMapping
%into the standard RAVEN model struct, applying the same defaults and
%post-processing the legacy parser does.

modelFields = legacyModelFields();
model = initModel(modelFields);

%Top-level scalar fields. The canonical schema lifts id/name to the top
%level; legacy callers sometimes wrote `version` at the top too. metaData
%(if present) carries provenance and defaultLB/defaultUB.
if isfield(doc,'id')
    model.id = stringify(doc.id);
end
if isfield(doc,'name')
    model.name = stringify(doc.name);
end
if isfield(doc,'version')
    model.version = stringify(doc.version);
end
if isfield(doc,'description')
    model.description = stringify(doc.description);
end

if isfield(doc,'metaData') && isstruct(doc.metaData)
    md = doc.metaData;
    %id/name in metaData take a back seat to top-level values (canonical)
    %but are still accepted (forward compat with legacy on-disk leakage).
    if isempty(model.id) && isfield(md,'id')
        model.id = stringify(md.id);
    end
    if isfield(md,'short_name') && isempty(model.id)
        model.id = stringify(md.short_name);
    end
    if isempty(model.name) && isfield(md,'name')
        model.name = stringify(md.name);
    end
    if isfield(md,'full_name') && isempty(model.description)
        model.description = stringify(md.full_name);
    end
    if isfield(md,'version') && isempty(model.version)
        model.version = stringify(md.version);
    end
    if isfield(md,'date')
        model.date = stringify(md.date);
    end
    if isfield(md,'taxonomy')
        model.annotation.taxonomy = stringify(md.taxonomy);
    end
    if isfield(md,'note')
        model.annotation.note = stringify(md.note);
    end
    if isfield(md,'description')
        model.annotation.note = stringify(md.description);
    end
    if isfield(md,'sourceUrl')
        model.annotation.sourceUrl = stringify(md.sourceUrl);
    end
    if isfield(md,'github') && ~isfield(model.annotation,'sourceUrl')
        model.annotation.sourceUrl = stringify(md.github);
    end
    if isfield(md,'givenName')
        model.annotation.givenName = stringify(md.givenName);
    end
    if isfield(md,'familyName')
        model.annotation.familyName = stringify(md.familyName);
    end
    if isfield(md,'authors')
        model.annotation.authorList = stringify(md.authors);
    end
    if isfield(md,'email')
        model.annotation.email = stringify(md.email);
    end
    if isfield(md,'organization')
        model.annotation.organization = stringify(md.organization);
    end
    if isfield(md,'defaultLB')
        model.annotation.defaultLB = str2double(stringify(md.defaultLB));
    end
    if isfield(md,'defaultUB')
        model.annotation.defaultUB = str2double(stringify(md.defaultUB));
    end
end

%Compartments (canonical: flat mapping `id -> name`).
if isfield(doc,'compartments') && isstruct(doc.compartments)
    cfields = fieldnames(doc.compartments);
    model.comps     = cell(numel(cfields),1);
    model.compNames = cell(numel(cfields),1);
    for k = 1:numel(cfields)
        model.comps{k,1}     = recoverKey(cfields{k});
        model.compNames{k,1} = stringify(doc.compartments.(cfields{k}));
    end
end

%Metabolites.
if isfield(doc,'metabolites') && iscell(doc.metabolites)
    nm = numel(doc.metabolites);
    model.mets        = cell(nm,1);
    model.metNames    = cell(nm,1);
    model.metComps    = cell(nm,1);
    model.metFormulas = cell(nm,1);
    model.metCharges  = cell(nm,1);
    model.inchis      = cell(nm,1);
    model.metSmiles   = cell(nm,1);
    model.metNotes    = cell(nm,1);
    model.metFrom     = cell(nm,1);
    model.metDeltaG   = cell(nm,1);
    model.metMiriams  = cell(nm,1);
    for k = 1:nm
        m = doc.metabolites{k};
        if ~isstruct(m), continue, end
        model.mets{k}        = getField(m,'id','');
        model.metNames{k}    = getField(m,'name','');
        model.metComps{k}    = getField(m,'compartment','');
        model.metFormulas{k} = getField(m,'formula','');
        model.metCharges{k}  = getField(m,'charge','');
        model.inchis{k}      = getField(m,'inchis','');
        if isfield(m,'annotation') && isstruct(m.annotation)
            ann = m.annotation;
            if isfield(ann,'smiles')
                model.metSmiles{k} = collapseList(ann.smiles);
                ann = rmfield(ann,'smiles');
            end
            if isfield(ann,'notes')
                model.metNotes{k} = collapseList(ann.notes);
                ann = rmfield(ann,'notes');
            end
            if isfield(ann,'metFrom')
                model.metFrom{k} = collapseList(ann.metFrom);
                ann = rmfield(ann,'metFrom');
            end
            if isfield(ann,'deltaG')
                model.metDeltaG{k} = collapseList(ann.deltaG);
                ann = rmfield(ann,'deltaG');
            end
            model.metMiriams{k} = annotationToMiriam(ann);
        end
        %Legacy callers placed these as top-level keys; honour both:
        if isempty(model.metSmiles{k}) && isfield(m,'smiles')
            model.metSmiles{k} = stringify(m.smiles);
        end
        if isempty(model.metNotes{k}) && isfield(m,'notes')
            model.metNotes{k} = stringify(m.notes);
        end
        if isempty(model.metFrom{k}) && isfield(m,'metFrom')
            model.metFrom{k} = stringify(m.metFrom);
        end
        if isempty(model.metDeltaG{k}) && isfield(m,'deltaG')
            model.metDeltaG{k} = stringify(m.deltaG);
        end
    end
end

%Reactions.
equations = {};
if isfield(doc,'reactions') && iscell(doc.reactions)
    nr = numel(doc.reactions);
    model.rxns                = cell(nr,1);
    model.rxnNames            = cell(nr,1);
    model.lb                  = cell(nr,1);
    model.ub                  = cell(nr,1);
    model.rev                 = cell(nr,1);
    model.grRules             = cell(nr,1);
    model.eccodes             = cell(nr,1);
    model.subSystems          = cell(nr,1);
    model.rxnReferences       = cell(nr,1);
    model.rxnNotes            = cell(nr,1);
    model.rxnFrom             = cell(nr,1);
    model.rxnConfidenceScores = cell(nr,1);
    model.rxnDeltaG           = cell(nr,1);
    model.rxnMiriams          = cell(nr,1);
    model.c                   = zeros(nr,1);
    for k = 1:nr
        r = doc.reactions{k};
        if ~isstruct(r), continue, end
        model.rxns{k}                = getField(r,'id','');
        model.rxnNames{k}            = getField(r,'name','');
        model.lb{k}                  = getField(r,'lower_bound','');
        model.ub{k}                  = getField(r,'upper_bound','');
        model.rev{k}                 = getField(r,'rev','');
        model.grRules{k}             = getField(r,'gene_reaction_rule','');
        model.rxnConfidenceScores{k} = getField(r,'confidence_score','');
        if isfield(r,'objective_coefficient')
            model.c(k) = str2double(stringify(r.objective_coefficient));
        end
        %Stoich (canonical: flat mapping; legacy: list of single-key maps
        %already flattened to a struct by parseMapping).
        if isfield(r,'metabolites') && isstruct(r.metabolites)
            stoichFields = fieldnames(r.metabolites);
            for s = 1:numel(stoichFields)
                metId = recoverKey(stoichFields{s});
                coeff = str2double(stringify(r.metabolites.(stoichFields{s})));
                equations(end+1,1:3) = {k, metId, coeff}; %#ok<AGROW>
            end
        end
        %Annotation block: ec-code, subsystem, references, notes, rxnFrom,
        %deltaG go to RAVEN's dedicated fields; everything else becomes
        %a miriam entry.
        if isfield(r,'annotation') && isstruct(r.annotation)
            ann = r.annotation;
            if isfield(ann,'ec_dash_code')
                model.eccodes{k} = strjoin(toCellList(ann.ec_dash_code),';');
                ann = rmfield(ann,'ec_dash_code');
            end
            if isfield(ann,'subsystem')
                model.subSystems{k} = toCellList(ann.subsystem);
                ann = rmfield(ann,'subsystem');
            end
            if isfield(ann,'references')
                model.rxnReferences{k} = collapseList(ann.references);
                ann = rmfield(ann,'references');
            end
            if isfield(ann,'notes')
                model.rxnNotes{k} = collapseList(ann.notes);
                ann = rmfield(ann,'notes');
            end
            if isfield(ann,'rxnFrom')
                model.rxnFrom{k} = collapseList(ann.rxnFrom);
                ann = rmfield(ann,'rxnFrom');
            end
            if isfield(ann,'deltaG')
                model.rxnDeltaG{k} = collapseList(ann.deltaG);
                ann = rmfield(ann,'deltaG');
            end
            model.rxnMiriams{k} = annotationToMiriam(ann);
        end
        %Top-level fields per the canonical schema. (Annotation rescue
        %above already covered the legacy-inside-annotation layout.)
        if isempty(model.eccodes{k}) && isfield(r,'eccodes')
            model.eccodes{k} = stringify(r.eccodes);
        end
        if isempty(model.rxnReferences{k}) && isfield(r,'references')
            model.rxnReferences{k} = stringify(r.references);
        end
        if (isempty(model.subSystems{k}) || (iscell(model.subSystems{k}) && isempty(model.subSystems{k}))) && isfield(r,'subsystem')
            sub = r.subsystem;
            if ischar(sub) || isstring(sub)
                model.subSystems{k} = {char(sub)};
            elseif iscell(sub)
                model.subSystems{k} = sub(:);
            end
        end
        if isempty(model.rxnNotes{k}) && isfield(r,'notes')
            model.rxnNotes{k} = stringify(r.notes);
        end
        if isempty(model.rxnFrom{k}) && isfield(r,'rxnFrom')
            model.rxnFrom{k} = stringify(r.rxnFrom);
        end
        if isempty(model.rxnDeltaG{k}) && isfield(r,'deltaG')
            model.rxnDeltaG{k} = stringify(r.deltaG);
        end
    end
end

%Genes.
if isfield(doc,'genes') && iscell(doc.genes)
    ng = numel(doc.genes);
    model.genes          = cell(ng,1);
    model.geneShortNames = cell(ng,1);
    model.proteins       = cell(ng,1);
    model.geneMiriams    = cell(ng,1);
    for k = 1:ng
        g = doc.genes{k};
        if ~isstruct(g), continue, end
        model.genes{k}          = getField(g,'id','');
        model.geneShortNames{k} = getField(g,'name','');
        model.proteins{k}       = getField(g,'protein','');
        if isfield(g,'annotation') && isstruct(g.annotation)
            model.geneMiriams{k} = annotationToMiriam(g.annotation);
        end
    end
end

%Gecko-light flag.
isGECKO = false;
if isfield(doc,'gecko_light')
    isGECKO = true;
    model.ec.geckoLight = isTruthy(doc.gecko_light);
end

%ec-rxns.
enzStoich = {};
if isfield(doc,'ec_dash_rxns') && iscell(doc.ec_dash_rxns)
    isGECKO = true;
    if ~isfield(model,'ec') || ~isfield(model.ec,'geckoLight')
        model.ec.geckoLight = false;
    end
    ne = numel(doc.ec_dash_rxns);
    model.ec.rxns    = cell(ne,1);
    model.ec.kcat    = cell(ne,1);
    model.ec.source  = cell(ne,1);
    model.ec.notes   = cell(ne,1);
    model.ec.eccodes = cell(ne,1);
    for k = 1:ne
        e = doc.ec_dash_rxns{k};
        if ~isstruct(e), continue, end
        model.ec.rxns{k}   = getField(e,'id','');
        model.ec.kcat{k}   = getField(e,'kcat','');
        model.ec.source{k} = getField(e,'source','');
        model.ec.notes{k}  = getField(e,'notes','');
        if isfield(e,'eccodes')
            model.ec.eccodes{k} = stringify(joinList(e.eccodes,';'));
        end
        if isfield(e,'enzymes') && isstruct(e.enzymes)
            sf = fieldnames(e.enzymes);
            for s = 1:numel(sf)
                enzId = recoverKey(sf{s});
                coeff = str2double(stringify(e.enzymes.(sf{s})));
                enzStoich(end+1,1:3) = {k, enzId, coeff}; %#ok<AGROW>
            end
        end
    end
end

%ec-enzymes.
if isfield(doc,'ec_dash_enzymes') && iscell(doc.ec_dash_enzymes)
    isGECKO = true;
    if ~isfield(model,'ec') || ~isfield(model.ec,'geckoLight')
        model.ec.geckoLight = false;
    end
    nEnz = numel(doc.ec_dash_enzymes);
    model.ec.genes    = cell(nEnz,1);
    model.ec.enzymes  = cell(nEnz,1);
    model.ec.mw       = cell(nEnz,1);
    model.ec.sequence = cell(nEnz,1);
    model.ec.concs    = cell(nEnz,1);
    for k = 1:nEnz
        e = doc.ec_dash_enzymes{k};
        if ~isstruct(e), continue, end
        model.ec.genes{k}    = getField(e,'genes','');
        model.ec.enzymes{k}  = getField(e,'enzymes','');
        model.ec.mw{k}       = getField(e,'mw','');
        model.ec.sequence{k} = getField(e,'sequence','');
        model.ec.concs{k}    = getField(e,'concs','');
    end
end

%Common finalisation: compartment lookup, stoich matrix, type coercions,
%empty-or-fill defaults, rxnGeneMat. Shared with parseLegacy.
model = finaliseModel(model, equations, enzStoich, isGECKO, verbose);
end

function out = getField(s,name,default)
if isfield(s,name)
    out = stringify(s.(name));
else
    out = default;
end
end

function out = stringify(v)
if iscell(v)
    if isempty(v)
        out = '';
    else
        out = char(v{1});
    end
elseif isstring(v)
    if numel(v) == 0
        out = '';
    else
        out = char(v(1));
    end
elseif isnumeric(v)
    if isnan(v)
        out = 'NaN';
    else
        out = sprintf('%.15g',v);
    end
elseif islogical(v)
    if v, out = 'true'; else, out = 'false'; end
else
    out = char(v);
end
end

function out = collapseList(v)
%Annotation values arrive as one-element lists in canonical YAML; collapse
%to a single string for RAVEN's scalar-string fields.
if iscell(v)
    parts = cellfun(@stringify,v,'UniformOutput',false);
    out = strjoin(parts,';');
elseif isstruct(v)
    %Shouldn't happen for these fields, but be permissive.
    out = '';
else
    out = stringify(v);
end
end

function out = toCellList(v)
if iscell(v)
    out = cellfun(@stringify,v,'UniformOutput',false);
elseif ischar(v) || isstring(v)
    out = {char(v)};
elseif isnumeric(v)
    out = {stringify(v)};
else
    out = {};
end
end

function out = joinList(v,sep)
out = strjoin(toCellList(v),sep);
end

function tf = isTruthy(v)
s = stringify(v);
tf = any(strcmpi(s,{'true','yes','1'}));
end

function key = recoverKey(safeKey)
%Reverse of sanitiseKey: strip the leading-digit guard, decode dots and
%dashes.
key = safeKey;
if startsWith(key,'kbeg_')
    key = key(6:end);
end
key = strrep(key,'_dot_','.');
key = strrep(key,'_dash_','-');
end

function miriam = annotationToMiriam(ann)
miriam = [];
if ~isstruct(ann)
    return
end
names  = {};
values = {};
fnames = fieldnames(ann);
for i = 1:numel(fnames)
    v = ann.(fnames{i});
    vals = toCellList(v);
    if isempty(vals)
        continue
    end
    displayName = recoverKey(fnames{i});
    for j = 1:numel(vals)
        names{end+1,1}  = displayName; %#ok<AGROW>
        values{end+1,1} = vals{j};     %#ok<AGROW>
    end
end
if isempty(names)
    miriam = [];
    return
end
miriam.name  = names;
miriam.value = values;
end

% --- Legacy-format parser (former monolithic body) ---------------------

function model = parseLegacy(line_raw, verbose)
% If entry is broken of multiple lines, concatenate. Assumes at least 6
% leading spaces to avoid metaData to be concatenated.
newLine=regexp(line_raw,'^ {6,}([\w\(\)].*)','tokenExtents');
brokenLine=find(~cellfun('isempty',newLine));
for i=flip(1:numel(brokenLine))
    extraLine = char(line_raw(brokenLine(i)));
    extraLine = extraLine(newLine{brokenLine(i)}{1}(1):end);
    line_raw{brokenLine(i)-1} = strjoin({line_raw{brokenLine(i)-1},extraLine},' ');
end
line_raw(brokenLine)=[];

line_key = regexprep(line_raw,'^ *-? ([^:]+)(:)($| .*)','$1');
line_key = regexprep(line_key,'(.*!!omap)|(---)|( {4,}.*)','');

line_value = regexprep(line_raw, '.*:$','');
line_value = regexprep(line_value, '[^":]+: "?(.+)"?$','$1');
line_value = regexprep(line_value, '(")|(^ {4,}- )','');
line_value(strcmp(line_value,'''''')) = {''};
line_value(strcmp(line_value,line_raw)) = {''};

model = initModel(legacyModelFields());

% If GECKO model
if any(contains(line_key,'geckoLight'))
    isGECKO=true;
    ecFields = {'geckoLight', false;...
                      'rxns', {};...
                      'kcat', {};...
                    'source', cell(0,0);...
                     'notes', cell(0,0);...
                   'eccodes', cell(0,0);...
                     'genes', cell(0,0);...
                   'enzymes', cell(0,0);...
                        'mw', cell(0,0);...
                  'sequence', cell(0,0);...
                     'concs', cell(0,0);...
                 'rxnEnzMat', []};
    for i=1:size(ecFields,1)
        model.ec.(ecFields{i,1})=ecFields{i,2};
    end
    ecGecko=cell(25000,2);      ecGeckoNo=1;
    enzStoich=cell(100000,3);   enzStoichNo=1;
else
    isGECKO=false;
end

section = 0;
metMiriams=cell(100000,3);   metMirNo=1;
rxnMiriams=cell(100000,3);   rxnMirNo=1;
geneMiriams=cell(100000,3);  genMirNo=1;
subSystems=cell(100000,2);   subSysNo=1;
eccodes=cell(100000,2);      ecCodeNo=1;
equations=cell(100000,3);    equatiNo=1;

for i=1:numel(line_key)
    tline_raw = line_raw{i};
    tline_key = line_key{i};
    tline_value = line_value{i};
    % import different sections
    switch tline_raw
        case '- metaData:'
            section = 1;
            if verbose
                fprintf('\t%d\n', section);
            end
            continue % Go to next line
        case '- metabolites:'
            section = 2;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- reactions:'
            section = 3;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- genes:'
            section = 4;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- compartments: !!omap'
            section = 5;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- ec-rxns:'
            section = 6;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
        case '- ec-enzymes:'
            section = 7;
            if verbose
                fprintf('\t%d\n', section);
            end
            pos=0;
            continue
    end

    % skip over empty keys
    if isempty(tline_raw) || (isempty(tline_key) && contains(tline_raw,'!!omap'))
        continue;
    end

    % import metaData
    if section == 1
        switch tline_key
            case {'short_name','id'} %short_name used by human-GEM
                model.id = tline_value;
            case 'name'
                model.name = tline_value;
            case 'full_name' %used by human-GEM
                model.description = tline_value;
            case 'version'
                model.version = tline_value;
            case 'date'
                model.date = tline_value;
            case 'taxonomy'
                model.annotation.taxonomy = tline_value;
            case {'description','note'} %description used by human-GEM
                model.annotation.note = tline_value;
            case 'github'
                model.annotation.sourceUrl = tline_value;
            case 'sourceUrl'
                model.annotation.sourceUrl = tline_value;
            case 'givenName'
                model.annotation.givenName = tline_value;
            case 'familyName'
                model.annotation.familyName = tline_value;
            case 'authors'
                model.annotation.authorList = tline_value;
            case 'email'
                model.annotation.email = tline_value;
            case 'organization'
                model.annotation.organization = tline_value;
            case 'geckoLight'
                if strcmp(tline_value,'true')
                    model.ec.geckoLight = true;
                end
        end; continue
    end

    % import metabolites:
    if section == 2
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'mets', tline_value,pos);
                readList=''; miriamKey='';
            case 'name'
                model = readFieldValue(model, 'metNames', tline_value, pos);
                readList=''; miriamKey='';
            case 'compartment'
                model = readFieldValue(model, 'metComps', tline_value, pos);
                readList=''; miriamKey='';
            case 'formula'
                model = readFieldValue(model, 'metFormulas', tline_value, pos);
                readList=''; miriamKey='';
            case 'charge'
                model = readFieldValue(model, 'metCharges', tline_value, pos);
                readList=''; miriamKey='';
            case 'notes'
                model = readFieldValue(model, 'metNotes', tline_value, pos);
                readList=''; miriamKey='';
            case 'inchis'
                model = readFieldValue(model, 'inchis', tline_value, pos);
                readList=''; miriamKey='';
            case 'smiles'
                model = readFieldValue(model, 'metSmiles', tline_value, pos);
                readList=''; miriamKey='';
            case 'deltaG'
                model = readFieldValue(model, 'metDeltaG', tline_value, pos);
                readList=''; miriamKey='';
            case 'metFrom'
                model = readFieldValue(model, 'metFrom', tline_value, pos);
                readList=''; miriamKey='';
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [metMiriams, miriamKey, metMirNo] = gatherAnnotation(pos,metMiriams,tline_key,tline_value,miriamKey,metMirNo);
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import reactions:
    if section == 3
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'rxns', tline_value,pos);
                readList=''; miriamKey='';
            case 'name'
                model = readFieldValue(model, 'rxnNames', tline_value, pos);
                readList=''; miriamKey='';
            case 'lower_bound'
                model.lb(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'upper_bound'
                model.ub(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'rev'
                model.rev(pos,1) = {tline_value};
                readList=''; miriamKey='';
            case 'gene_reaction_rule'
                model = readFieldValue(model, 'grRules', tline_value, pos);
                readList=''; miriamKey='';
            case 'rxnNotes'
                model = readFieldValue(model, 'rxnNotes', tline_value, pos);
                readList=''; miriamKey='';
            case 'rxnFrom'
                model = readFieldValue(model, 'rxnFrom', tline_value, pos);
                readList=''; miriamKey='';
            case 'deltaG'
                model = readFieldValue(model, 'rxnDeltaG', tline_value, pos);
                readList=''; miriamKey='';
            case 'objective_coefficient'
                model.c(pos,1) = 1;
                readList=''; miriamKey='';
            case 'references'
                model = readFieldValue(model, 'rxnReferences', tline_value, pos);
                readList=''; miriamKey='';
            case 'confidence_score'
                model = readFieldValue(model, 'rxnConfidenceScores', tline_value, pos);
                readList=''; miriamKey='';
            case 'eccodes'
                if isempty(tline_value)
                    readList = 'eccodes';
                else
                    eccodes(ecCodeNo,1:2)={pos,tline_value};
                    ecCodeNo=ecCodeNo+1;
                end
            case 'subsystem'
                if isempty(tline_value)
                    readList = 'subsystem';
                else
                    subSystems(subSysNo,1:2)={pos,tline_value};
                    subSysNo=subSysNo+1;
                end
            case 'metabolites'
                readList = 'equation';
            case 'annotation'
                readList = 'annotation';

            otherwise
                switch readList
                    case 'eccodes'
                        eccodes(ecCodeNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        ecCodeNo=ecCodeNo+1;
                    case 'subsystem'
                        subSystems(subSysNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        subSysNo=subSysNo+1;
                    case 'annotation'
                        [rxnMiriams, miriamKey,rxnMirNo] = gatherAnnotation(pos,rxnMiriams,tline_key,tline_value,miriamKey,rxnMirNo);
                    case 'equation'
                        coeff = sscanf(tline_value,'%f');
                        equations(equatiNo,1:3)={pos,tline_key,coeff};
                        equatiNo=equatiNo+1;
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import genes:
    if section == 4
        switch tline_key
            case 'id'
                pos = pos + 1;
                model = readFieldValue(model, 'genes', tline_value, pos);
                readList = '';
                miriamKey = '';
            case 'name'
                model = readFieldValue(model, 'geneShortNames', tline_value, pos);
            case 'protein'
                model = readFieldValue(model, 'proteins', tline_value, pos);
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [geneMiriams, miriamKey,genMirNo] = gatherAnnotation(pos,geneMiriams,tline_key,tline_value,miriamKey,genMirNo);
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import compartments:
    if section == 5
        model.comps(end+1,1) = {tline_key};
        model.compNames(end+1,1) = {tline_value};
    end

    % import ec reaction info
    if section == 6
        switch tline_key
            case 'id'
                pos = pos + 1;
                model.ec = readFieldValue(model.ec, 'rxns', tline_value, pos);
                readList='';
            case 'kcat'
                model.ec = readFieldValue(model.ec, 'kcat', tline_value, pos);
                readList='';
            case 'source'
                model.ec = readFieldValue(model.ec, 'source', tline_value, pos);
                readList='';
            case 'notes'
                model.ec = readFieldValue(model.ec, 'notes', tline_value, pos);
                readList='';
            case 'eccodes'
                if isempty(tline_value)
                    readList = 'eccodes';
                else
                    ecGecko(ecGeckoNo,1:2)={pos,tline_value};
                    ecGeckoNo=ecGeckoNo+1;
                end
            case 'enzymes'
                readList = 'enzStoich';
            otherwise
                switch readList
                    case 'eccodes'
                        ecGecko(ecGeckoNo,1:2)={pos,regexprep(tline_value,'^ +- "?(.*)"?$','$1')};
                        ecGeckoNo=ecGeckoNo+1;
                    case 'enzStoich'
                        coeff = sscanf(tline_value,'%f');
                        enzStoich(enzStoichNo,1:3)={pos,tline_key,coeff};
                        enzStoichNo=enzStoichNo+1;
                    otherwise
                        error(['Unknown entry in yaml file: ' tline_raw])
                end
        end; continue
    end

    % import ec enzyme info
    if section == 7
        switch tline_key
            case 'genes'
                pos = pos + 1;
                model.ec = readFieldValue(model.ec, 'genes', tline_value, pos);
            case 'enzymes'
                model.ec = readFieldValue(model.ec, 'enzymes', tline_value, pos);
            case 'mw'
                model.ec = readFieldValue(model.ec, 'mw', tline_value, pos);
            case 'sequence'
                model.ec = readFieldValue(model.ec, 'sequence', tline_value, pos);
            case 'concs'
                model.ec = readFieldValue(model.ec, 'concs', tline_value, pos);
            otherwise
                error(['Unknown entry in yaml file: ' tline_raw])
        end; continue
    end
end

%Parse annotations
if ~isempty(metMiriams)
    locs = cell2mat(metMiriams(:,1));
    for i=unique(locs)'
        model.metMiriams{i,1}.name=metMiriams(locs==i,2);
        model.metMiriams{i,1}.value=metMiriams(locs==i,3);
    end
end
if ~isempty(rxnMiriams)
    locs = cell2mat(rxnMiriams(:,1));
    for i=unique(locs)'
        model.rxnMiriams{i,1}.name=rxnMiriams(locs==i,2);
        model.rxnMiriams{i,1}.value=rxnMiriams(locs==i,3);
    end
end
if ~isempty(geneMiriams)
    locs = cell2mat(geneMiriams(:,1));
    for i=unique(locs)'
        model.geneMiriams{i,1}.name=geneMiriams(locs==i,2);
        model.geneMiriams{i,1}.value=geneMiriams(locs==i,3);
    end
end

%Parse subSystems
if ~isempty(subSystems)
    locs = cell2mat(subSystems(:,1));
    for i=unique(locs)'
        model.subSystems{i,1}=subSystems(locs==i,2);
    end
end

%Parse ec-codes
if ~isempty(eccodes)
    locs = cell2mat(eccodes(:,1));
    for i=unique(locs)'
        eccodesCat=strjoin(eccodes(locs==i,2),';');
        model.eccodes{i,1}=eccodesCat;
    end
    emptyEc=cellfun('isempty',model.eccodes);
    model.eccodes(emptyEc)={''};
end

% Pack ecGecko (eccodes for ec-rxns) so the shared finaliser can consume it.
ecGeckoOut = {};
if isGECKO
    ecGeckoOut = ecGecko;
    enzStoichOut = enzStoich;
else
    enzStoichOut = {};
end

model = finaliseModel(model, equations, enzStoichOut, isGECKO, verbose, ecGeckoOut);
end

% --- Shared finalisation -----------------------------------------------

function modelFields = legacyModelFields()
modelFields =   {'id',char();...
               'name',char();...
        'description',char();...
            'version',char();...
               'date',char();...
         'annotation',struct();...
               'rxns',{};...
           'rxnNames',{};...
               'mets',{};...
           'metNames',{};...
                  'S',sparse([]);...
                 'lb',{};... %Changed to double in the end.
                 'ub',{};... %Changed to double in the end.
                'rev',{};... %Changed to double in the end.
                  'c',[];...
                  'b',cell(0,0);... %Changed to double in the end.
              'genes',cell(0,0);...
            'grRules',cell(0,0);...
         'rxnGeneMat',sparse([]);...
           'rxnComps',cell(0,0);... %Changed to double in the end.
         'subSystems',cell(0,0);...
            'eccodes',cell(0,0);...
         'rxnMiriams',cell(0,0);...
          'rxnDeltaG',{};... %Changed to double in the end.
           'rxnNotes',cell(0,0);...
      'rxnReferences',cell(0,0);...
'rxnConfidenceScores',cell(0,0);...
           'metComps',cell(0,0);... %Changed to double in the end.
             'inchis',cell(0,0);...
          'metSmiles',cell(0,0);...
        'metFormulas',cell(0,0);...
         'metMiriams',cell(0,0);...
          'metDeltaG',{};... %Changed to double in the end.
         'metCharges',cell(0,0);... %Changed to double in the end.
           'metNotes',cell(0,0);...
              'comps',cell(0,0);...
          'compNames',cell(0,0);...
        'compOutside',cell(0,0);...
          'geneComps',cell(0,0);... %Changed to double in the end.
        'geneMiriams',cell(0,0);...
     'geneShortNames',cell(0,0);...
       'proteins',cell(0,0);...
      'unconstrained',cell(0,0);... %Changed to double in the end.
            'metFrom',cell(0,0);...
            'rxnFrom',cell(0,0)};
end

function model = initModel(modelFields)
model = struct();
for i=1:size(modelFields,1)
    model.(modelFields{i,1})=modelFields{i,2};
end
end

function model = finaliseModel(model, equations, enzStoich, isGECKO, verbose, ecGecko)
if nargin < 6
    ecGecko = {};
end

if verbose
    fprintf('\nimporting completed\nfollow-up processing...');
end
[~, model.metComps] = ismember(model.metComps, model.comps);
[~, model.geneComps] = ismember(model.geneComps, model.comps);
[~, model.rxnComps] = ismember(model.rxnComps, model.comps);

% Fill S-matrix
if ~isempty(equations)
    rxnIdx = cellfun('isempty', equations(:,1));
    equations(rxnIdx,:) = '';
    if ~isempty(equations)
        rxnIdx = cell2mat(equations(:,1));
        [~,metIdx] = ismember(equations(:,2),model.mets);
        coeffs = cell2mat(equations(:,3));
        model.S=sparse(max(metIdx),max(rxnIdx));
        linearIndices = sub2ind([max(metIdx), max(rxnIdx)],metIdx,rxnIdx);
        model.S(linearIndices) = coeffs;
    end
end

% Convert strings to numeric
model.metCharges = str2double(model.metCharges);
model.lb = str2double(model.lb);
model.ub = str2double(model.ub);
model.rxnConfidenceScores = str2double(model.rxnConfidenceScores);
model.b = zeros(length(model.mets),1);
model.metDeltaG = str2double(model.metDeltaG);
model.rxnDeltaG = str2double(model.rxnDeltaG);

% Fill some other fields
model.annotation.defaultLB = min(model.lb);
model.annotation.defaultUB = max(model.ub);
if numel(model.lb)<numel(model.rxns) %No LB reported = min
    model.lb(end+1:numel(model.rxns)-numel(model.lb),1) = double(model.annotation.defaultLB);
end
if numel(model.ub)<numel(model.rxns) %No UB reported = max
    model.ub(end+1:numel(model.rxns)-numel(model.ub),1) = double(model.annotation.defaultUB);
end
if ~all(cellfun('isempty',model.rev))
    model.rev = str2double(model.rev);
else
    model.rev = [];
end
if numel(model.rev)<numel(model.rxns) %No rev reported, assume from LB and UB
    model.rev(end+1:numel(model.rxns)-numel(model.rev),1) = double(model.lb<0 & model.ub>0);
end

% Remove empty fields, otherwise fill to correct length
% Reactions
for i={'rxnNames','grRules','eccodes','rxnNotes','rxnReferences',...
       'rxnFrom','subSystems','rxnMiriams'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'rxns');
end
for i={'c'} % Zeros
   model = emptyOrFill(model,i{1},0,'rxns',true);
end
for i={'rxnConfidenceScores','rxnDeltaG'} % NaNs
   model = emptyOrFill(model,i{1},NaN,'rxns');
end
for i={'rxnComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'rxns');
end
% Metabolites
for i={'metNames','inchis','metFormulas','metMiriams','metFrom','metSmiles','metNotes'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'mets');
end
for i={'metCharges','unconstrained'} % Zeros
   model = emptyOrFill(model,i{1},0,'mets');
end
for i={'metDeltaG'} % NaNs
    model = emptyOrFill(model,i{1},NaN,'mets');
 end
for i={'metComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'mets');
end
% Genes
for i={'geneMiriams','geneShortNames','proteins'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'genes');
end
for i={'geneComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'genes');
end
% Comps
for i={'compNames'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'comps');
end
for i={'compOutside'} % First comp
   model = emptyOrFill(model,i{1},model.comps{1},'comps');
end

if isfield(model,'subSystems')
% If all entries are 1x1, then flatten
    if all(cellfun(@(x) numel(x) <= 1, model.subSystems))
        model.subSystems = transpose([model.subSystems{:}]);
    end
end

% Make rxnGeneMat fields and map to the existing model.genes field
if isfield(model,'grRules')
    [genes, rxnGeneMat] = getGenesFromGrRules(model.grRules);
    model.rxnGeneMat = sparse(numel(model.rxns),numel(model.genes));
    [~,geneOrder] = ismember(genes,model.genes);
    if any(geneOrder == 0)
        error(['The grRules includes the following gene(s), that are not in '...
            'the list of model genes: ', genes{~geneOrder}])
    end
    model.rxnGeneMat(:,geneOrder) = rxnGeneMat;
end

% Finalize GECKO model
if isGECKO
    % Fill in empty fields and empty entries
    for i={'kcat','source','notes','eccodes'} % Even keep empty
        model.ec = emptyOrFill(model.ec,i{1},{''},'rxns',true);
    end
    for i={'enzymes','mw','sequence'}
        model.ec = emptyOrFill(model.ec,i{1},{''},'genes',true);
    end
    model.ec = emptyOrFill(model.ec,'concs',{'NaN'},'genes',true);
    model.ec = emptyOrFill(model.ec,'kcat',{'0'},'genes',true);
    % Change string to double
    for i={'kcat','mw','concs'}
        if isfield(model.ec,i{1})
            model.ec.(i{1}) = str2double(model.ec.(i{1}));
        end
    end
    % Fill rxnEnzMat
    if ~isempty(enzStoich)
        rxnIdx              = cellfun('isempty', enzStoich(:,1));
        enzStoich(rxnIdx,:) = '';
        if ~isempty(enzStoich)
            rxnIdx              = cell2mat(enzStoich(:,1));
            [~,enzIdx]          = ismember(enzStoich(:,2),model.ec.enzymes);
            coeffs              = cell2mat(enzStoich(:,3));
            model.ec.rxnEnzMat  = zeros(numel(model.ec.rxns), numel(model.ec.genes));
            linearIndices       = sub2ind([numel(model.ec.rxns), numel(model.ec.genes)], rxnIdx, enzIdx);
            model.ec.rxnEnzMat(linearIndices) = coeffs;
        end
    end
    %Parse ec-codes (legacy multi-line eccodes path).
    if ~isempty(ecGecko)
        locs = cell2mat(ecGecko(:,1));
        for i=unique(locs)'
            ecGeckoCat=strjoin(ecGecko(locs==i,2),';');
            model.ec.eccodes{i,1}=ecGeckoCat;
        end
        emptyEc=cellfun('isempty',model.ec.eccodes);
        model.ec.eccodes(emptyEc)={''};
    end
end

if verbose
    fprintf(' Done!\n');
end
end

function model = emptyOrFill(model,field,emptyEntry,type,keepEmpty)
if nargin<5
    keepEmpty=false;
end
if isnumeric(emptyEntry)
    emptyCells=isempty(model.(field));
else
    emptyCells=cellfun('isempty',model.(field));
end
if all(emptyCells) && ~keepEmpty
    model = rmfield(model, field);
elseif numel(model.(field))<numel(model.(type))
    model.(field)(end+1:numel(model.(type)),1)=emptyEntry;
end
end

function model = readFieldValue(model, fieldName, value, pos)
if numel(model.(fieldName))<pos-1
    model.(fieldName)(end+1:pos,1) = {''};
end
model.(fieldName)(pos,1) = {value};
end

function [miriams,miriamKey,entryNumber] = gatherAnnotation(pos,miriams,key,value,miriamKey,entryNumber)
if isempty(key)
    key=miriamKey;
else
    miriamKey=key;
end
if ~isempty(value)
    miriams(entryNumber,1:3) = {pos, key, strip(value)};
    entryNumber = entryNumber + 1;
end
end
