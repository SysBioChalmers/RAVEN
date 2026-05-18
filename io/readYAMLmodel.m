function model=readYAMLmodel(fileName, verbose)
% readYAMLmodel
%   Reads a model from a YAML file.
%
%   Two YAML layouts are recognised and handled automatically:
%
%   1. The cobrapy-style layout that writeYAMLmodel now produces.
%      This is what cobrapy and any tools built on top of it (e.g.
%      Escher, Memote) read and write.
%
%   2. The older RAVEN layout (files starting with `---` followed by
%      `!!omap`). These are still accepted so that YAML files written
%      by earlier RAVEN/GECKO releases keep loading without manual
%      conversion.
%
%   In both cases the function returns a regular RAVEN model struct.
%
%   Input:
%   fileName    path to the YAML file. If empty, a dialog window opens
%               so the user can pick the file.
%   verbose     if true, print progress messages while reading
%               (optional, default false).
%
%   Output:
%   model       a RAVEN model structure.
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

% Pick which parser to use. The older RAVEN files start with `---`
% and use `!!omap` tags; the newer cobrapy-style files do not.
if isLegacyFormat(line_raw)
    model = parseYAMLLegacy(line_raw, verbose);
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
%Parser for the cobrapy-style YAML layout. The file is read into a
%nested struct/cell tree first, then translated into a standard RAVEN
%model struct.

if verbose
    fprintf('Parsing YAML (cobrapy layout)...\n');
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
            %`-`. Capture that first key (either as an inline scalar, or
            %as a nested block when the value comes on the following
            %lines), then run the continuation loop to pick up any
            %remaining keys of the same element.
            element = struct();
            safeKey = sanitiseKey(key);
            if hasInline
                element.(safeKey) = parseScalar(value);
                idx = idx + 1;
            else
                nextIdx = idx + 1;
                if nextIdx > numel(lines)
                    element.(safeKey) = '';
                    idx = nextIdx;
                else
                    [nIndent, nContent] = splitIndent(lines{nextIdx});
                    if nIndent <= baseIndent + 1
                        %No nested block; key has an empty value.
                        element.(safeKey) = '';
                        idx = nextIdx;
                    else
                        if startsWith(nContent,'-')
                            [grand, idx] = parseList(lines, nextIdx, nIndent);
                        else
                            [grand, idx] = parseMapping(lines, nextIdx, nIndent);
                        end
                        element.(safeKey) = grand;
                    end
                end
            end
            %Continuation: subsequent lines indented at (baseIndent + 2)
            %belong to this same element (its remaining keys).
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
                    %Reached the next list element belonging to a parent
                    %list (e.g. annotation values). Stop here.
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
if s(1) ~= '"' && s(1) ~= ''''
    commentStart = regexp(s,'\s#','once');
    if ~isempty(commentStart)
        s = strtrim(s(1:commentStart));
    end
end
%Strip surrounding quotes (either ".." or '..'). YAML uses both, and
%cobrapy/geckopy emit single quotes for empty strings via PyYAML.
if numel(s) >= 2 && ((s(1) == '"' && s(end) == '"') || (s(1) == '''' && s(end) == ''''))
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

modelFields = legacyYAMLmodelFields();
model = initYAMLmodel(modelFields);

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

%Metabolites. String-cell fields are seeded with '' so unset entries
%end up char-empty (matching the legacy pipeline), not [] empty-double.
%Numeric-cell fields (metCharges, metDeltaG) stay as cell(nm,1) so that
%after str2double their unset entries become NaN, which finaliseYAMLmodel
%treats as "all empty" and drops the field entirely.
if isfield(doc,'metabolites') && iscell(doc.metabolites)
    nm = numel(doc.metabolites);
    model.mets        = repmat({''},nm,1);
    model.metNames    = repmat({''},nm,1);
    model.metComps    = repmat({''},nm,1);
    model.metFormulas = repmat({''},nm,1);
    model.metCharges  = cell(nm,1);
    model.inchis      = repmat({''},nm,1);
    model.metSmiles   = repmat({''},nm,1);
    model.metNotes    = repmat({''},nm,1);
    model.metFrom     = repmat({''},nm,1);
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
    model.rxns                = repmat({''},nr,1);
    model.rxnNames            = repmat({''},nr,1);
    model.lb                  = cell(nr,1);
    model.ub                  = cell(nr,1);
    model.rev                 = cell(nr,1);
    model.grRules             = repmat({''},nr,1);
    model.eccodes             = repmat({''},nr,1);
    model.subSystems          = cell(nr,1);
    model.rxnReferences       = repmat({''},nr,1);
    model.rxnNotes            = repmat({''},nr,1);
    model.rxnFrom             = repmat({''},nr,1);
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
    model.genes          = repmat({''},ng,1);
    model.geneShortNames = repmat({''},ng,1);
    model.proteins       = repmat({''},ng,1);
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
    %String fields pre-seeded with '' so any unset entries match the
    %char-empty convention the legacy pipeline produces (cell(N,1)
    %would have left them as the [] empty-double default).
    model.ec.rxns    = repmat({''},ne,1);
    model.ec.kcat    = cell(ne,1);
    model.ec.source  = repmat({''},ne,1);
    model.ec.notes   = repmat({''},ne,1);
    model.ec.eccodes = repmat({''},ne,1);
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
    model.ec.genes    = repmat({''},nEnz,1);
    model.ec.enzymes  = repmat({''},nEnz,1);
    model.ec.mw       = cell(nEnz,1);
    model.ec.sequence = repmat({''},nEnz,1);
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
model = finaliseYAMLmodel(model, equations, enzStoich, isGECKO, verbose);
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

