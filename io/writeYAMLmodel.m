function writeYAMLmodel(model,fileName,preserveQuotes,sortIds)
% writeYAMLmodel
%   Writes a model to a YAML file.
%
%   The format follows the same layout used by cobrapy (the Python
%   COBRA toolbox), so the resulting file can also be opened by
%   cobrapy and any tool built on top of it (e.g. Escher, Memote).
%   Enzyme-constrained extras (the `ec-rxns`, `ec-enzymes`,
%   `gecko_light` and `metaData` sections) sit alongside the regular
%   model content; cobrapy ignores them, while readYAMLmodel reads
%   them back in full.
%
%   Older RAVEN YAML files (the ones with `---` / `!!omap` headers)
%   are no longer produced, but the companion reader readYAMLmodel
%   still loads them, so existing files keep working.
%
%   model           a model structure
%   fileName        name that the file will have. A dialog window will
%                   open if no file name is specified.
%   preserveQuotes  if true, string values are written with surrounding
%                   double quotes (logical, default=true)
%   sortIds         if true, metabolites, reactions, genes and
%                   compartments are written in alphabetical order;
%                   otherwise they are kept in their original order
%                   (logical, default=false)
%
% Usage: writeYAMLmodel(model,fileName,preserveQuotes,sortIds)
if nargin<2|| isempty(fileName)
    [fileName, pathName] = uiputfile({'*.yml;*.yaml'}, 'Select file for model export',[model.id '.yml']);
    if fileName == 0
        error('You should provide a file location')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);

if nargin < 3
    preserveQuotes = true;
end
if nargin < 4
    sortIds = false;
end
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

%Pre-compute simplified miriam annotation tables (namespace -> values):
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

%Model id and name go at the top level, alongside cobrapy's other
%top-level keys:
if isfield(model,'id') && ~isempty(model.id)
    fprintf(fid,'id: %s\n',quoteIfNeeded(model.id,preserveQuotes));
else
    fprintf(fid,'id: %s\n',quoteIfNeeded('blankID',preserveQuotes));
end
if isfield(model,'name') && ~isempty(model.name)
    fprintf(fid,'name: %s\n',quoteIfNeeded(model.name,preserveQuotes));
end
if isfield(model,'version') && ~isempty(model.version)
    fprintf(fid,'version: %s\n',quoteIfNeeded(model.version,preserveQuotes));
end

%Optional metaData block (provenance such as date, author, taxonomy;
%cobrapy ignores it but it is preserved on round-trip):
writeMetadata(model,fid,preserveQuotes);

%Top-level gecko_light flag (native boolean):
if isfield(model,'ec')
    if isfield(model.ec,'geckoLight') && model.ec.geckoLight
        fprintf(fid,'gecko_light: true\n');
    else
        fprintf(fid,'gecko_light: false\n');
    end
end

%Compartments (flat mapping):
if isfield(model,'comps') && ~isempty(model.comps)
    fprintf(fid,'compartments:\n');
    for i = 1:length(model.comps)
        compName = '';
        if isfield(model,'compNames') && i <= numel(model.compNames) && ~isempty(model.compNames{i})
            compName = model.compNames{i};
        end
        fprintf(fid,'  %s: %s\n',model.comps{i},quoteIfNeeded(compName,preserveQuotes));
    end
end

%Metabolites:
fprintf(fid,'metabolites:\n');
for i = 1:length(model.mets)
    writeListItemField(model, fid, 'mets',        'txt', i, 'id',          preserveQuotes, true)
    writeListItemField(model, fid, 'metNames',    'txt', i, 'name',        preserveQuotes, false)
    writeListItemField(model, fid, 'metComps',    'txt', i, 'compartment', preserveQuotes, false)
    writeListItemField(model, fid, 'metFormulas', 'txt', i, 'formula',     preserveQuotes, false)
    writeListItemField(model, fid, 'metCharges',  'num', i, 'charge',      preserveQuotes, false)
    writeListItemField(model, fid, 'inchis',      'txt', i, 'inchis',      preserveQuotes, false)
    writeListItemField(model, fid, 'metDeltaG',   'num', i, 'deltaG',      preserveQuotes, false)
    writeListItemField(model, fid, 'metNotes',    'txt', i, 'notes',       preserveQuotes, false)
    writeListItemField(model, fid, 'metFrom',     'txt', i, 'metFrom',     preserveQuotes, false)
    %smiles + miriam namespace fields land inside the annotation block:
    writeAnnotation(model, fid, 'metabolite', i, preserveQuotes)
end

%Reactions:
fprintf(fid,'reactions:\n');
for i = 1:length(model.rxns)
    writeListItemField(model, fid, 'rxns',                 'txt', i, 'id',                    preserveQuotes, true)
    writeListItemField(model, fid, 'rxnNames',             'txt', i, 'name',                  preserveQuotes, false)
    writeStoichMapping(model, fid, i)
    writeListItemField(model, fid, 'lb',                   'num', i, 'lower_bound',           preserveQuotes, false)
    writeListItemField(model, fid, 'ub',                   'num', i, 'upper_bound',           preserveQuotes, false)
    writeListItemField(model, fid, 'grRules',              'txt', i, 'gene_reaction_rule',    preserveQuotes, false)
    if isfield(model,'c') && i <= numel(model.c) && model.c(i) ~= 0
        writeListItemField(model, fid, 'c',                'num', i, 'objective_coefficient', preserveQuotes, false)
    end
    writeSubsystem(model, fid, i, preserveQuotes)
    writeListItemField(model, fid, 'rxnConfidenceScores',  'num', i, 'confidence_score',      preserveQuotes, false)
    writeListItemField(model, fid, 'rxnReferences',        'txt', i, 'references',            preserveQuotes, false)
    writeListItemField(model, fid, 'rxnFrom',              'txt', i, 'rxnFrom',               preserveQuotes, false)
    writeListItemField(model, fid, 'rxnNotes',             'txt', i, 'notes',                 preserveQuotes, false)
    writeListItemField(model, fid, 'rxnDeltaG',            'num', i, 'deltaG',                preserveQuotes, false)
    %ec-code + miriam namespace fields land inside the annotation block:
    writeAnnotation(model, fid, 'reaction', i, preserveQuotes)
end

%Genes:
if isfield(model,'genes')
    fprintf(fid,'genes:\n');
    for i = 1:length(model.genes)
        writeListItemField(model, fid, 'genes',          'txt', i, 'id',      preserveQuotes, true)
        writeListItemField(model, fid, 'geneShortNames', 'txt', i, 'name',    preserveQuotes, false)
        writeListItemField(model, fid, 'proteins',       'txt', i, 'protein', preserveQuotes, false)
        writeAnnotation(model, fid, 'gene', i, preserveQuotes)
    end
end

%EC-model:
if isfield(model,'ec')
    fprintf(fid,'ec-rxns:\n');
    for i = 1:length(model.ec.rxns)
        writeListItemField(model.ec, fid, 'rxns',    'txt', i, 'id',      preserveQuotes, true)
        writeListItemField(model.ec, fid, 'kcat',    'num', i, 'kcat',    preserveQuotes, false)
        writeListItemField(model.ec, fid, 'source',  'txt', i, 'source',  preserveQuotes, false)
        writeListItemField(model.ec, fid, 'notes',   'txt', i, 'notes',   preserveQuotes, false)
        writeEcEccodes(model.ec, fid, i, preserveQuotes)
        writeEnzymeMapping(model.ec, fid, i)
    end

    fprintf(fid,'ec-enzymes:\n');
    for i = 1:length(model.ec.genes)
        writeListItemField(model.ec, fid, 'genes',    'txt', i, 'genes',    preserveQuotes, true)
        writeListItemField(model.ec, fid, 'enzymes',  'txt', i, 'enzymes',  preserveQuotes, false)
        writeListItemField(model.ec, fid, 'mw',       'num', i, 'mw',       preserveQuotes, false)
        writeListItemField(model.ec, fid, 'sequence', 'txt', i, 'sequence', preserveQuotes, false)
        writeListItemField(model.ec, fid, 'concs',    'num', i, 'concs',    preserveQuotes, false)
    end
end

%Close file:
fclose(fid);

end

% --- Helpers -----------------------------------------------------------

function s = quoteIfNeeded(value,preserveQuotes)
%Returns a YAML scalar representation of a string. Empty -> '""'.
if isempty(value)
    s = '""';
    return
end
if ~ischar(value) && ~isstring(value)
    value = char(value);
end
if preserveQuotes
    s = ['"' char(value) '"'];
else
    s = char(value);
end
end

function writeListItemField(model,fid,fieldName,type,pos,name,preserveQuotes,isFirst)
%Write one "field: value" line for the current list item. When isFirst
%is true the leading "- " bullet of the item is emitted as well.
%Missing or empty values are silently skipped; the caller only sets
%isFirst=true for the id field, which is always present.
if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps stores the compartment as a numeric index into
        %model.comps; the YAML file uses the compartment id string.
        if pos > numel(model.metComps) || model.metComps(pos) == 0
            value = '';
        else
            value = model.comps{model.metComps(pos)};
        end
        if isempty(value)
            return
        end
        emitLine(fid,name,quoteIfNeeded(value,preserveQuotes),isFirst);
        return
    end

    field = model.(fieldName);
    if strcmp(type,'txt')
        if pos > numel(field)
            return
        end
        value = field{pos};
        if isempty(value)
            return
        end
        emitLine(fid,name,quoteIfNeeded(value,preserveQuotes),isFirst);
    elseif strcmp(type,'num')
        if pos > numel(field)
            return
        end
        v = full(field(pos));
        if isnan(v)
            %For mw/concs we write NaN explicitly as `.nan` so the value
            %survives the round trip (NaN is their natural "no value"
            %marker). For other numeric fields NaN just means "no
            %value": skip the line and let the reader fall back to its
            %default.
            if ismember(fieldName,{'mw','concs'})
                emitLine(fid,name,'.nan',isFirst);
            end
            return
        end
        emitLine(fid,name,sprintf('%.15g',v),isFirst);
    end
end
end

function emitLine(fid,name,value,isFirst)
if isFirst
    fprintf(fid,'- %s: %s\n',name,value);
else
    fprintf(fid,'  %s: %s\n',name,value);
end
end

function writeStoichMapping(model,fid,rxnPos)
%Emit a reaction's metabolites as a flat mapping:
%    metabolites:
%      met_id_1: coeff
%      met_id_2: coeff
if ~isfield(model,'S')
    return
end
col = model.S(:,rxnPos);
mask = col ~= 0;
if ~any(mask)
    return
end
mets   = model.mets(mask);
coeffs = full(col(mask));
[mets,order] = sort(mets);
coeffs       = coeffs(order);
fprintf(fid,'  metabolites:\n');
for i = 1:numel(mets)
    fprintf(fid,'    %s: %.15g\n',mets{i},coeffs(i));
end
end

function writeEnzymeMapping(ec,fid,rxnPos)
%Emit an ec-rxn's enzymes as a flat mapping:
%    enzymes:
%      enzyme_id_1: stoich
%      enzyme_id_2: stoich
if ~isfield(ec,'rxnEnzMat') || isempty(ec.rxnEnzMat)
    return
end
row = ec.rxnEnzMat(rxnPos,:);
mask = row ~= 0;
if ~any(mask)
    return
end
enzymes = ec.enzymes(mask);
coeffs  = row(mask);
[enzymes,order] = sort(enzymes);
coeffs          = coeffs(order);
fprintf(fid,'  enzymes:\n');
for i = 1:numel(enzymes)
    fprintf(fid,'    %s: %.15g\n',enzymes{i},coeffs(i));
end
end

function writeEcEccodes(ec,fid,pos,preserveQuotes)
%ec-rxns[].eccodes: write as a single quoted string when there is exactly
%one code, otherwise as a list.
if ~isfield(ec,'eccodes') || pos > numel(ec.eccodes) || isempty(ec.eccodes{pos})
    return
end
codes = strsplit(strrep(ec.eccodes{pos},' ',''),';');
codes = strip(codes);
codes = codes(~cellfun('isempty',codes));
if isempty(codes)
    return
end
if numel(codes) == 1
    fprintf(fid,'  eccodes: %s\n',quoteIfNeeded(codes{1},preserveQuotes));
else
    fprintf(fid,'  eccodes:\n');
    for i = 1:numel(codes)
        fprintf(fid,'  - %s\n',quoteIfNeeded(codes{i},preserveQuotes));
    end
end
end

function writeAnnotation(model,fid,kind,pos,preserveQuotes)
%Write the `annotation:` block for a metabolite, reaction or gene.
%It only contains identifier annotations (MIRIAM cross-references,
%plus SMILES for metabolites and EC numbers for reactions); other
%per-entry fields (notes, references, subsystem, ...) are written
%as top-level keys on the entry itself.
miriam = struct('names',{{}},'values',{{}});
switch kind
    case 'metabolite'
        if isfield(model,'newMetMiriams') && size(model.newMetMiriams,1) >= pos
            miriam = harvestMiriam(model.newMetMiriams,model.newMetMiriamNames,pos);
        end
    case 'reaction'
        if isfield(model,'newRxnMiriams') && size(model.newRxnMiriams,1) >= pos
            miriam = harvestMiriam(model.newRxnMiriams,model.newRxnMiriamNames,pos);
        end
    case 'gene'
        if isfield(model,'newGeneMiriams') && size(model.newGeneMiriams,1) >= pos
            miriam = harvestMiriam(model.newGeneMiriams,model.newGeneMiriamNames,pos);
        end
end

extras = {};
switch kind
    case 'metabolite'
        if isfield(model,'metSmiles') && pos <= numel(model.metSmiles) && ~isempty(model.metSmiles{pos})
            extras(end+1,:) = {'smiles',{model.metSmiles{pos}}};
        end
    case 'reaction'
        if isfield(model,'eccodes') && pos <= numel(model.eccodes) && ~isempty(model.eccodes{pos})
            codes = strsplit(strrep(model.eccodes{pos},' ',''),';');
            codes = strip(codes);
            codes = codes(~cellfun('isempty',codes));
            if ~isempty(codes)
                extras(end+1,:) = {'ec-code',codes};
            end
        end
end

if isempty(miriam.names) && isempty(extras)
    return
end

fprintf(fid,'  annotation:\n');
%Miriam entries first (sorted by namespace for reproducibility):
[sortedNames,order] = sort(miriam.names);
for i = 1:numel(sortedNames)
    values = miriam.values{order(i)};
    writeAnnotationValues(fid,sortedNames{i},values,preserveQuotes);
end
%Then the extras, in the order added above:
for i = 1:size(extras,1)
    writeAnnotationValues(fid,extras{i,1},extras{i,2},preserveQuotes);
end
end

function writeSubsystem(model,fid,pos,preserveQuotes)
%Emit a reaction's subsystem as a top-level key. A single-element cell
%collapses to a quoted scalar; multi-element cells go out as a YAML list.
if ~isfield(model,'subSystems') || pos > numel(model.subSystems)
    return
end
sub = model.subSystems{pos};
if isempty(sub)
    return
end
if ischar(sub) || isstring(sub)
    sub = {char(sub)};
end
if ~iscell(sub)
    return
end
sub = sub(~cellfun('isempty',sub));
if isempty(sub)
    return
end
if numel(sub) == 1
    fprintf(fid,'  subsystem: %s\n',quoteIfNeeded(sub{1},preserveQuotes));
else
    fprintf(fid,'  subsystem:\n');
    for i = 1:numel(sub)
        fprintf(fid,'  - %s\n',quoteIfNeeded(sub{i},preserveQuotes));
    end
end
end

function out = harvestMiriam(matrix,names,pos)
out = struct('names',{{}},'values',{{}});
for j = 1:size(matrix,2)
    cell_ij = matrix{pos,j};
    if isempty(cell_ij)
        continue
    end
    vals = strsplit(cell_ij,'; ');
    vals = strip(vals);
    vals = vals(~cellfun('isempty',vals));
    if isempty(vals)
        continue
    end
    out.names{end+1,1}  = names{j};
    out.values{end+1,1} = vals;
end
end

function writeAnnotationValues(fid,key,values,preserveQuotes)
%Write each annotation value on its own line, so the resulting block
%is a YAML list of strings (which is how cobrapy expects annotations).
fprintf(fid,'    %s:\n',key);
for i = 1:numel(values)
    fprintf(fid,'    - %s\n',quoteIfNeeded(values{i},preserveQuotes));
end
end

function writeMetadata(model,fid,preserveQuotes)
%Write the optional `metaData:` block. It collects provenance fields
%(date, taxonomy, author, ...); the model id, name and gecko_light
%flag live as top-level keys, not inside this block.
fprintf(fid,'metaData:\n');
fprintf(fid, '  date: %s\n', quoteIfNeeded(datestr(now,29),preserveQuotes));  % 29=YYYY-MM-DD
if isfield(model,'annotation')
    a = model.annotation;
    fields = {'defaultLB','defaultUB','givenName','familyName', ...
              'authors','email','organization','taxonomy','note','sourceUrl'};
    for k = 1:numel(fields)
        f = fields{k};
        if isfield(a,f) && ~isempty(a.(f))
            v = a.(f);
            if ischar(v) || isstring(v)
                fprintf(fid,'  %s: %s\n',f,quoteIfNeeded(char(v),preserveQuotes));
            else
                fprintf(fid,'  %s: %g\n',f,v);
            end
        end
    end
end
end
