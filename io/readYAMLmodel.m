function model=readYAMLmodel(yamlFilename, verbose)
% readYAMLmodel
%   Reads a yaml file matching (roughly) the cobrapy yaml structure
%
%   Input:
%   yamlFile    a model file in yaml file format
%   verbose     set as true to monitor progress (opt, default false)
%
%   Output:
%   model       a model structure
%
%   Usage: model = readYAMLmodel(yamlFilename, verbose)

if nargin < 2
    verbose = false;
end

if ~(exist(yamlFilename,'file')==2)
    error('Yaml file %s cannot be found', string(yamlFilename));
end

if verLessThan('matlab','9.9') %readlines introduced 2020b
    fid=fopen(yamlFilename);
    line_raw=cell(1000000,1);
    while ~feof(fid)
        line_raw{i}=fgetl(fid);
        i=i+1;
    end
    line_raw(i:end)=[];
    line_raw=string(line_raw);
else
    line_raw=readlines(yamlFilename);
end

% If entry is broken of multiple lines, concatenate. Assumes at least 6
% leading spaces to avoid metaData to be concatenated.
newLine=regexp(line_raw,'^ {6,}([\w\(\)].*)','tokens');
brokenLine=find(~cellfun(@isempty,newLine));
for i=1:numel(brokenLine)
    line_raw{brokenLine(i)-1} = strjoin({line_raw{brokenLine(i)-1},char(newLine{brokenLine(i)}{1})},' ');
end
line_raw(brokenLine)=[];

line_key = regexprep(line_raw,'^ *-? ([^:]+)(:).*','$1');
line_key = regexprep(line_key,'(.*!!omap)|(---)|( {4,}.*)','');

line_value = regexprep(line_raw, '.*:$','');
line_value = regexprep(line_value, '[^":]+: "?(.+)"?$','$1');
line_value = regexprep(line_value, '(")|(^ {4,}- )','');

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
           'rxnNotes',cell(0,0);...
      'rxnReferences',cell(0,0);...
'rxnConfidenceScores',cell(0,0);...
           'metComps',cell(0,0);... %Changed to double in the end.
             'inchis',cell(0,0);...
          'metSmiles',cell(0,0);...
        'metFormulas',cell(0,0);...
         'metMiriams',cell(0,0);...
         'metCharges',cell(0,0);... %Changed to double in the end.
              'comps',cell(0,0);...
          'compNames',cell(0,0);...
        'compOutside',cell(0,0);...
          'geneComps',cell(0,0);... %Changed to double in the end.
        'geneMiriams',cell(0,0);...
     'geneShortNames',cell(0,0);...
      'unconstrained',cell(0,0);... %Changed to double in the end.
            'metFrom',cell(0,0);...
            'rxnFrom',cell(0,0)};
for i=1:size(modelFields,1)
    model.(modelFields{i,1})=modelFields{i,2};
end

section = 0;
metMiriams=cell(25000,3);   metMirNo=1;
rxnMiriams=cell(25000,3);   rxnMirNo=1;
geneMiriams=cell(25000,3);  genMirNo=1;
subSystems=cell(25000,2);   subSysNo=1;
eccodes=cell(25000,2);      ecCodeNo=1;
equations=cell(100000,3);   equatiNo=1;
% metMiriams=cell(0,3);   metMirNo=1;
% rxnMiriams=cell(0,3);   rxnMirNo=1;
% geneMiriams=cell(0,3);  genMirNo=1;
% subSystems=cell(0,2);   subSysNo=1;
% eccodes=cell(0,2);      ecCodeNo=1;

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
            case 'inchis'
                model = readFieldValue(model, 'inchis', tline_value, pos);
                readList=''; miriamKey='';
            case 'smiles'
                model = readFieldValue(model, 'metSmiles', tline_value, pos);
                readList=''; miriamKey='';                
            case 'metFrom'
                model = readFieldValue(model, 'metFrom', tline_value, pos);
                readList=''; miriamKey='';
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [metMiriams, miriamKey] = gatherAnnotation(pos,metMiriams,tline_key,tline_value,miriamKey,metMirNo);
                        metMirNo = metMirNo + 1;
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
                        rxnMirNo=rxnMirNo+1;
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
            case 'annotation'
                readList = 'annotation';
            otherwise
                switch readList
                    case 'annotation'
                        [geneMiriams, miriamKey] = gatherAnnotation(pos,geneMiriams,tline_key,tline_value,miriamKey,genMirNo);
                        genMirNo = genMirNo + 1;
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
    emptyEc=cellfun(@isempty,model.eccodes);
    model.eccodes(emptyEc)={''};
end

% follow-up data processing
if verbose
    fprintf('\nimporting completed\nfollow-up processing...');
end
[~, model.metComps] = ismember(model.metComps, model.comps);
[~, model.geneComps] = ismember(model.geneComps, model.comps);
[~, model.rxnComps] = ismember(model.rxnComps, model.comps);

% Fill S-matrix
rxnIdx = cellfun(@isempty, equations(:,1));
equations(rxnIdx,:) = '';
rxnIdx = cell2mat(equations(:,1));
[~,metIdx] = ismember(equations(:,2),model.mets);
coeffs = cell2mat(equations(:,3));
model.S=sparse(max(metIdx),max(rxnIdx));
linearIndices = sub2ind([max(metIdx), max(rxnIdx)],metIdx,rxnIdx);
model.S(linearIndices) = coeffs;

% Convert strings to numeric
model.metCharges = str2double(model.metCharges);
model.lb = str2double(model.lb);
model.ub = str2double(model.ub);
model.rxnConfidenceScores = str2double(model.rxnConfidenceScores);
model.b = zeros(length(model.mets),1);

% Fill some other fields
model.annotation.defaultLB = min(model.lb);
model.annotation.defaultUB = max(model.ub);
if numel(model.lb)<numel(model.rxns) %No LB reported = min
    model.lb(end+1:numel(model.rxns)-numel(model.lb),1) = double(model.annotation.defaultLB);
end
if numel(model.ub)<numel(model.rxns) %No UB reported = max
    model.ub(end+1:numel(model.rxns)-numel(model.ub),1) = double(model.annotation.defaultUB);
end
if ~all(cellfun(@isempty,model.rev))
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
   model = emptyOrFill(model,i{1},0,'rxns');
end
for i={'rxnConfidenceScores'} % NaNs
   model = emptyOrFill(model,i{1},NaN,'rxns');
end
for i={'rxnComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'rxns');
end
% Metabolites
for i={'metNames','inchis','metFormulas','metMiriams','metFrom','metSmiles'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'mets');
end
for i={'metCharges','unconstrained'} % Zeros
   model = emptyOrFill(model,i{1},0,'mets');
end
for i={'metComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'mets');
end
% Genes
for i={'geneMiriams','geneShortNames'} % Empty strings
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
% Single fields are kept, even if empty
% for i={'description','name','version','date','annotation'}
%     if isempty(model.(i{1}))
%         model = rmfield(model,i{1});
%     end
% end

% Make rxnGeneMat fields
[genes, rxnGeneMat] = getGenesFromGrRules(model.grRules, model.genes);
if isequal(sort(genes), sort(model.genes))
    model.rxnGeneMat = rxnGeneMat;
    model.genes = genes;
else
    error('The gene list and grRules are inconsistent.');
end

if verbose
    fprintf(' Done!\n');
end
end

function model = emptyOrFill(model,field,emptyEntry,type)
if isnumeric(emptyEntry)
    emptyCells=isempty(model.(field));
else
    emptyCells=cellfun(@isempty,model.(field));
end
if all(emptyCells)
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

function [miriams, miriamKey,entryNumber] = gatherAnnotation(pos,miriams,key,value,miriamKey,entryNumber)
if isempty(key)
    key=miriamKey;
else
    miriamKey=key;
end
if ~isempty(value)
    miriams(entryNumber,1:3) = {pos, key, strip(value)};
else
    entryNumber = entryNumber - 1;
end
end
