function model = parseYAMLLegacy(line_raw, verbose)
% parseYAMLLegacy
%   Read a RAVEN model from the older YAML layout (outer `---` /
%   `!!omap` wrapper, `id`/`name`/`geckoLight` nested under `metaData`,
%   per-entry `!!omap` tags, scalar annotation values).
%
%   This is the original line-by-line regex parser that shipped with
%   RAVEN before the switch to the cobrapy-style YAML layout. It is
%   called from readYAMLmodel when isLegacyFormat detects the older
%   shape; new files take the canonical path in readYAMLmodel itself.
%
% Input:
%   line_raw    cell/string array of YAML lines, as read by readlines
%               (or fgetl on older MATLAB releases).
%   verbose     true to print which section is being parsed.
%
% Output:
%   model       a RAVEN model structure (same shape as the canonical
%               parser produces).

% If an entry was broken across multiple lines, concatenate. Assumes
% at least 6 leading spaces so metaData lines are not touched.
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

model = initYAMLmodel(legacyYAMLmodelFields());

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

model = finaliseYAMLmodel(model, equations, enzStoichOut, isGECKO, verbose, ecGeckoOut);
end

% --- Legacy-only helpers ----------------------------------------------

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
