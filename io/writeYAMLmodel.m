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
% preserveQuotes : logical
%     if all string values should be wrapped in double quotes. cobrapy
%     emits quotes only where YAML requires them, so the default is false
%     (matches cobrapy / raven-python) (default false).
% sortIds : logical
%     if metabolites, reactions, genes and compartments should be sorted
%     alphabetically by their identifier, otherwise they are kept in their
%     original order (default false).
%
% Examples
% --------
%     writeYAMLmodel(model,fileName,preserveQuotes,sortIds);
p=parseRAVENargs(varargin, {'fileName',[]; 'preserveQuotes',false; 'sortIds',false});
fileName=p.fileName; preserveQuotes=p.preserveQuotes; sortIds=p.sortIds;
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
writeMetadata(model,fid,preserveQuotes);

%Metabolites:
% Field order matches cobrapy + raven_python.io.yaml:
%   id, name, compartment, charge, formula, notes, annotation,
%   then RAVEN-only extras (inchis, deltaG, metFrom).
% SMILES goes inside the annotation block (cobrapy convention), not at
% metabolite top level — the reader still accepts top-level `smiles:`
% for backward compatibility with older yeast-GEM files.
fprintf(fid,'- metabolites:\n');
for i = 1:length(model.mets)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'mets',        'txt', i, '  - id',          preserveQuotes)
    writeField(model, fid, 'metNames',    'txt', i, '  - name',        preserveQuotes)
    writeField(model, fid, 'metComps',    'txt', i, '  - compartment', preserveQuotes)
    writeField(model, fid, 'metCharges',  'num', i, '  - charge',      preserveQuotes)
    writeField(model, fid, 'metFormulas', 'txt', i, '  - formula',     preserveQuotes)
    writeField(model, fid, 'metNotes',    'txt', i, '  - notes',       preserveQuotes)
    writeAnnotation(model, fid, 'met',    i,      preserveQuotes)
    writeField(model, fid, 'inchis',      'txt', i, '  - inchis',      preserveQuotes)
    writeField(model, fid, 'metDeltaG',   'num', i, '  - deltaG',      preserveQuotes)
    writeField(model, fid, 'metFrom',     'txt', i, '  - metFrom',     preserveQuotes)
end

%Reactions:
% Field order matches cobrapy + raven_python.io.yaml:
%   id, name, metabolites, lower_bound, upper_bound, gene_reaction_rule,
%   objective_coefficient, subsystem, notes, annotation (which carries
%   the EC numbers under `ec-code`, the cobrapy/geckopy convention),
%   then RAVEN-only extras (references, rxnFrom, deltaG,
%   confidence_score). The notes key is the canonical `notes` (no
%   longer `rxnNotes`); the reader still accepts the legacy key.
fprintf(fid,'- reactions:\n');
for i = 1:length(model.rxns)
    fprintf(fid,'    - !!omap\n');
    writeField(model, fid, 'rxns',                 'txt', i, '  - id',                    preserveQuotes)
    writeField(model, fid, 'rxnNames',             'txt', i, '  - name',                  preserveQuotes)
    writeField(model, fid, 'S',                    'txt', i, '  - metabolites',           preserveQuotes)
    writeField(model, fid, 'lb',                   'num', i, '  - lower_bound',           preserveQuotes)
    writeField(model, fid, 'ub',                   'num', i, '  - upper_bound',           preserveQuotes)
    writeField(model, fid, 'grRules',              'txt', i, '  - gene_reaction_rule',    preserveQuotes)
    if model.c(i)~=0
        writeField(model, fid, 'c',                'num', i, '  - objective_coefficient', preserveQuotes)
    end
    writeField(model, fid, 'subSystems',           'txt', i, '  - subsystem',             preserveQuotes)
    writeField(model, fid, 'rxnNotes',             'txt', i, '  - notes',                 preserveQuotes)
    writeAnnotation(model, fid, 'rxn',             i,                                     preserveQuotes)
    writeField(model, fid, 'rxnReferences',        'txt', i, '  - references',            preserveQuotes)
    writeField(model, fid, 'rxnFrom',              'txt', i, '  - rxnFrom',               preserveQuotes)
    writeField(model, fid, 'rxnDeltaG',            'num', i, '  - deltaG',                preserveQuotes)
    writeField(model, fid, 'rxnConfidenceScores',  'num', i, '  - confidence_score',      preserveQuotes)
end

%Genes:
if isfield(model,'genes')
    fprintf(fid,'- genes:\n');
    for i = 1:length(model.genes)
        fprintf(fid,'    - !!omap\n');
        writeField(model, fid, 'genes',          'txt', i, '  - id',         preserveQuotes)
        writeField(model, fid, 'geneShortNames', 'txt', i, '  - name',       preserveQuotes)
        writeField(model, fid, 'proteins',   'txt', i, '  - protein',    preserveQuotes)
        writeField(model, fid, 'geneMiriams',    'txt', i, '  - annotation', preserveQuotes)
    end
end

%Compartments:
fprintf(fid,'- compartments: !!omap\n');
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames',   'txt', i, ['- ' model.comps{i}], preserveQuotes)
    writeField(model, fid, 'compMiriams', 'txt', i, '- annotation',             preserveQuotes)
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
        writeField(model.ec, fid, 'rxns',      'txt', i, '- id',      preserveQuotes)
        writeField(model.ec, fid, 'kcat',      'num', i, '- kcat',    preserveQuotes)
        writeField(model.ec, fid, 'source',    'txt', i, '- source',  preserveQuotes)
        writeField(model.ec, fid, 'notes',     'txt', i, '- notes',   preserveQuotes)
        writeField(model.ec, fid, 'eccodes',   'txt', i, '- eccodes', preserveQuotes)
        writeField(model.ec, fid, 'rxnEnzMat', 'txt', i, '- enzymes', preserveQuotes)
    end

    fprintf(fid,'- ec-enzymes:\n');
    for i = 1:length(model.ec.genes)
        fprintf(fid,'  - !!omap\n');
        writeField(model.ec, fid, 'genes',    'txt', i, '- genes',    preserveQuotes)
        writeField(model.ec, fid, 'enzymes',  'txt', i, '- enzymes',  preserveQuotes)
        writeField(model.ec, fid, 'mw',       'num', i, '- mw',       preserveQuotes)
        writeField(model.ec, fid, 'sequence', 'txt', i, '- sequence', preserveQuotes)
        writeField(model.ec, fid, 'concs',    'num', i, '- concs',    preserveQuotes)
    end
end

%Close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name,preserveQuotes)
%Writes a new line in the yaml file if the field exists and the field is
%not empty at the correspoinding position. It's recursive for some fields
%(metMiriams, rxnMiriams, and S)

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end
    
    field = model.(fieldName);
    
    if strcmp(fieldName,'metMiriams')
        if ~isempty(model.metMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newMetMiriams,2)
                %'i' represents the different miriam names, e.g.
                %kegg.compound or chebi
                if ~isempty(model.newMetMiriams{pos,i})
                    %As during the following writeField call the value of
                    %'i' would be lost, it is temporarily concatenated to
                    %'name' parameter, which will be edited later
                    writeField(model, fid, 'newMetMiriams', 'txt', pos, ['      - ' model.newMetMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'rxnMiriams')
        if ~isempty(model.rxnMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newRxnMiriams,2)
                if ~isempty(model.newRxnMiriams{pos,i})
                    writeField(model, fid, 'newRxnMiriams', 'txt', pos, ['      - ' model.newRxnMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'geneMiriams')
        if ~isempty(model.geneMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newGeneMiriams,2)
                if ~isempty(model.newGeneMiriams{pos,i})
                    writeField(model, fid, 'newGeneMiriams', 'txt', pos, ['      - ' model.newGeneMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'compMiriams')
        if ~isempty(model.compMiriams{pos})
            fprintf(fid,'    %s: !!omap\n',name);
            for i=1:size(model.newCompMiriams,2)
                if ~isempty(model.newCompMiriams{pos,i})
                    writeField(model, fid, 'newCompMiriams', 'txt', pos, ['      - ' model.newCompMiriamNames{i} '_' sprintf('%d',i)], preserveQuotes)
                end
            end
        end
        
    elseif strcmp(fieldName,'S')
        %S: create header & write each metabolite in a new line. Reactions
        %with no metabolites emit `metabolites: !!omap []` (the flow-style
        %empty omap cobrapy uses) so the file remains a valid YAML 1.2
        %document.
        if sum(field(:,pos) ~= 0) > 0
            fprintf(fid,'    %s: !!omap\n',name);
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            %Sort metabolites:
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'num', i, ['      - ' model.mets{i}], preserveQuotes)
            end
        else
            fprintf(fid,'    %s: !!omap []\n',name);
        end

    elseif strcmp(fieldName,'rxnEnzMat')
        %S: create header & write each enzyme in a new line
        fprintf(fid,'    %s: !!omap\n',name);
        if sum(field(pos,:) ~= 0) > 0
            model.enzymes = model.enzymes(field(pos,:) ~= 0);
            model.coeffs  = field(pos,field(pos,:) ~= 0);
            %Sort metabolites:
            [model.enzymes,order] = sort(model.enzymes);
            model.coeffs          = model.coeffs(order);
            for i = 1:length(model.enzymes)
                writeField(model, fid, 'coeffs',  'num', i, ['    - ' model.enzymes{i}], preserveQuotes)
            end
        end        

    elseif sum(strcmp({'subSystems','newMetMiriams','newRxnMiriams','newGeneMiriams','newCompMiriams','eccodes'},fieldName)) > 0
        %eccodes/rxnNotes: if 1 write in 1 line, if more create header and list
        if strcmp(fieldName,'subSystems')
            % The reader collapses an all-singleton subSystems field to
            % a char column; defend against that (and against length
            % mismatches caused by partial subsystem coverage) so the
            % writer doesn't crash on shorter-than-rxns subsystem lists.
            if iscell(field)
                if pos > numel(field); return; end
                list = field{pos};
            else
                if pos > size(field, 1); return; end
                list = field(pos, :);
            end
            if isempty(list)
                return
            end
        elseif strcmp(fieldName,'newMetMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newMetMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newRxnMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newRxnMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newGeneMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newGeneMiriams{pos,index},'; ');
        elseif strcmp(fieldName,'newCompMiriams')
            index = str2double(regexprep(name,'^.+_',''));
            name  = regexprep(name,'_\d+$','');
            list = strsplit(model.newCompMiriams{pos,index},'; ');
        elseif ~isempty(field{pos})
            list = strrep(field{pos},' ','');
            list = strsplit(list,';');
        else
            return % empty, needs no line in file
        end
        list=strip(list);

        if length(list) == 1 && ~strcmp(list{1},'') && ~strcmp(fieldName,'subSystems')
            if preserveQuotes
                list = ['"' list{1} '"'];
            end
            if iscell(list)
                list=list{1};
            end
            fprintf(fid,'    %s: %s\n',name,list);
        elseif ischar(list) && strcmp(fieldName,'subSystems')
            if preserveQuotes
                list = ['"' list '"'];
            end
            fprintf(fid,'    %s: %s\n',name,list);            
        elseif length(list) > 1 || strcmp(fieldName,'subSystems')
            if preserveQuotes
                for j=1:numel(list)
                    list{j} = ['"' list{j} '"'];
                end
            end
            fprintf(fid,'    %s:\n',name);
            for i = 1:length(list)
                fprintf(fid,'%s        - %s\n',regexprep(name,'(^\s*).*','$1'),list{i});
            end
        end
        
    elseif sum(pos) > 0
        %All other fields:
        if strcmp(type,'txt')
            value = field{pos};
            if ~isempty(value)
                if preserveQuotes || needsYamlQuoting(value)
                    value = ['"', escapeForDoubleQuoted(value), '"'];
                end
            end
        elseif strcmp(type,'num')
            if isnan(field(pos))
                value = [];
            else
                value = formatNumber(full(field(pos)));
            end
        end
        if ~isempty(value)
            fprintf(fid,'    %s: %s\n',name,value);
        end
    end
end
end

function writeMetadata(model, fid, preserveQuotes)
% Writes the metaData block. Honors preserveQuotes so the rest of the
% file (which defaults to no surrounding quotes for cobra parity) stays
% consistent. The `date` field is preserved when the model carries one
% (model.date), so round-trips don't churn on every write; if absent
% it's filled with the current date.

fprintf(fid, '- metaData: !!omap\n');
emitMetaField(fid, 'id',           valueOrDefault(model,'id','blankID'),    preserveQuotes);
emitMetaField(fid, 'name',         valueOrDefault(model,'name','blankName'),preserveQuotes);
if isfield(model,'version')
    emitMetaField(fid, 'version', model.version, preserveQuotes);
end
if isfield(model,'date') && ~isempty(model.date)
    dateValue = model.date;
else
    dateValue = datestr(now, 29); %#ok<DATST>  % 29 = yyyy-mm-dd
end
emitMetaField(fid, 'date', dateValue, preserveQuotes);
if isfield(model,'annotation')
    annoFields = {'defaultLB','defaultUB','givenName','familyName', ...
                  'authors','email','organization','taxonomy','note','sourceUrl'};
    for k = 1:numel(annoFields)
        f = annoFields{k};
        if isfield(model.annotation, f)
            emitMetaField(fid, f, model.annotation.(f), preserveQuotes);
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

function emitMetaField(fid, key, value, preserveQuotes)
% Emit one `  - key: value` line inside the metaData omap block.
if islogical(value)
    if value; value = 'true'; else; value = 'false'; end
end
if isnumeric(value)
    value = formatNumber(double(value));
elseif ~ischar(value) && ~isstring(value)
    value = char(value);
end
if preserveQuotes
    fprintf(fid, '  - %s: "%s"\n', key, value);
else
    fprintf(fid, '  - %s: %s\n', key, value);
end
end

function tf = needsYamlQuoting(s)
% A defensive subset of when a plain (unquoted) YAML scalar would be
% misparsed: leading YAML indicator chars, or any character that
% triggers flow-style collection / mapping parsing. Matches what
% ruamel.yaml's "round-trip" emitter quotes automatically.
if isempty(s)
    tf = false;
    return;
end
% Leading-character cases that turn into a flow indicator / tag / etc.
first = s(1);
if any(first == '[]{},&*!|>%@`#')
    tf = true; return;
end
if first == '-' || first == '?' || first == ':'
    tf = true; return;
end
% In-string cases: ': ' (key/value confusion), ' #' (comment), any flow
% bracket, leading or trailing whitespace, or anything outside ASCII
% printable.
if contains(s, ': ') || contains(s, ' #') || any(ismember(s, '[]{},'))
    tf = true; return;
end
if ~strcmp(strip(string(s)), string(s))
    tf = true; return;
end
% YAML reserves certain tokens (true, false, null, …) as bare scalars
% reading as booleans / nulls. Quote when the entire value matches.
if any(strcmpi(s, {'true','false','null','yes','no','on','off','~'}))
    tf = true; return;
end
tf = false;
end

function out = escapeForDoubleQuoted(s)
% Escape backslashes and double quotes for emission inside YAML's
% double-quoted style. Conservative — full YAML escapes (Unicode etc.)
% are out of scope for the strings model curators normally use.
out = strrep(s, '\', '\\');
out = strrep(out, '"', '\"');
end

function s = formatNumber(x)
% Format a finite number the way cobrapy / Python's float repr would —
% so whole-number floats round-trip as "1000.0", not "1000". Matches the
% ruamel.yaml output used by raven_python.io.yaml.write_yaml_model.
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

function writeAnnotation(model, fid, kind, pos, preserveQuotes)
% Emit the per-entry `annotation` block, fusing MIRIAM cross-references
% with the non-MIRIAM cobrapy-style annotation keys: `smiles` for
% metabolites and `ec-code` (EC numbers) for reactions. cobrapy and
% geckopy read both from inside `annotation` (not as top-level entry
% keys), so this helper keeps the YAML aligned. `ec-code` is emitted as a
% list — matching cobrapy / raven-python and geckopy, which read it as
% list[str] — even when there is a single code.
switch kind
    case 'met'
        miriamsField     = 'metMiriams';
        miriamNamesField = 'newMetMiriamNames';
        miriamValsField  = 'newMetMiriams';
        extraName        = 'smiles';
        extraField       = 'metSmiles';
        extraIsList      = false;
    case 'rxn'
        miriamsField     = 'rxnMiriams';
        miriamNamesField = 'newRxnMiriamNames';
        miriamValsField  = 'newRxnMiriams';
        extraName        = 'ec-code';
        extraField       = 'eccodes';
        extraIsList      = true;
    otherwise
        error('writeAnnotation:unsupportedKind', 'Unsupported kind: %s', kind);
end

hasMiriams = isfield(model, miriamsField) && ~isempty(model.(miriamsField){pos});
hasExtra   = isfield(model, extraField) && pos <= numel(model.(extraField)) ...
             && ~isempty(model.(extraField){pos});

if ~hasMiriams && ~hasExtra
    return;
end

fprintf(fid, '      - annotation: !!omap\n');
if hasMiriams
    % Flat fprintf over the extractMiriam intermediate (the block header
    % is already written above).
    miriamNames  = model.(miriamNamesField);
    miriamValues = model.(miriamValsField);
    for j = 1:size(miriamValues, 2)
        v = miriamValues{pos, j};
        if isempty(v); continue; end
        list = strip(strsplit(strrep(v, ' ', ''), ';'));
        if numel(list) == 1
            fprintf(fid, '          - %s: %s\n', miriamNames{j}, ...
                quoteIfNeeded(list{1}, preserveQuotes));
        else
            fprintf(fid, '          - %s:\n', miriamNames{j});
            for k = 1:numel(list)
                fprintf(fid, '              - %s\n', ...
                    quoteIfNeeded(list{k}, preserveQuotes));
            end
        end
    end
end
if hasExtra
    if extraIsList
        % EC numbers: a ;-joined string -> a block list under `ec-code`
        % (always a list, matching cobrapy/raven-python/geckopy).
        codes = strip(strsplit(strrep(model.(extraField){pos}, ' ', ''), ';'));
        codes = codes(~cellfun('isempty', codes));
        if ~isempty(codes)
            fprintf(fid, '          - %s:\n', extraName);
            for k = 1:numel(codes)
                fprintf(fid, '              - %s\n', ...
                    quoteIfNeeded(codes{k}, preserveQuotes));
            end
        end
    else
        fprintf(fid, '          - %s: %s\n', extraName, ...
            quoteIfNeeded(model.(extraField){pos}, preserveQuotes));
    end
end
end

function out = quoteIfNeeded(value, preserveQuotes)
% Quote a YAML scalar exactly when the surrounding writer would have
% lost it as a flow sequence / boolean / null otherwise. Conservative
% wrapper around needsYamlQuoting that respects preserveQuotes=true.
if preserveQuotes || needsYamlQuoting(value)
    out = ['"', escapeForDoubleQuoted(value), '"'];
else
    out = value;
end
end
