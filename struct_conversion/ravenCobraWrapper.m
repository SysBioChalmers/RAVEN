function newModel=ravenCobraWrapper(model)
% ravenCobraWrapper
%   Converts between RAVEN and COBRA structures
%
%   Input: model          a RAVEN/COBRA-compatible model structure
%
%   Ouput: newModel       a COBRA/RAVEN-compatible model structure
%   
%   This function is a bidirectional tool to convert between RAVEN and
%   COBRA structures. It recognises COBRA structure by checking field
%   'rules' existense, which is only found in COBRA Toolbox structure.
%
%   NOTE: During RAVEN -> COBRA -> RAVEN conversion cycle the following
%   fields are lost: annotation, compOutside, compMiriams, rxnComps,
%   geneComps, unconstrained. Boundary metabolites are lost, because COBRA
%   structure does not involve boundary metabolites, so they are removed
%   using simplifyModel before RAVEN -> COBRA conversion. The field 'rev'
%   is also partially lost, but during COBRA -> RAVEN conversion it's
%   reconstructed based on lower bound reaction values
%
%   NOTE: During COBRA -> RAVEN -> COBRA conversion cycle the following
%   fields are lost: geneEntrezID, metSmiles, modelVersion,
%   proteinNames, proteins
%
%   NOTE: The information about mandatory RAVEN fields was taken from
%   checkModelStruct function, whereas the corresponding information about
%   COBRA fields was fetched from verifyModel function
%
%   Usage: newModel=ravenCobraWrapper(model)

if isfield(model,'rules')
    isRaven=false;
else
    isRaven=true;
end

[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

% Load COBRA field information
fid             = fopen([ravenPath filesep 'struct_conversion' filesep 'COBRA_structure_fields.csv']); % Taken from https://github.com/opencobra/cobratoolbox/blob/develop/src/base/io/definitions/COBRA_structure_fields.csv
fieldFile       = textscan(fid,repmat('%s',1,15),'Delimiter','\t','HeaderLines',1);
dbFields        = ~cellfun(@isempty,fieldFile{5}); % Only keep fields with database annotations that should be translated to xxxMiriams
dbFields        = dbFields & ~contains(fieldFile{1},{'metInChIString','metKEGGID','metPubChemID','rxnECNumbers'});
COBRAnamespace  = fieldFile{5}(dbFields);
COBRAnamespace  = regexprep(COBRAnamespace,';.*',''); % Only keep first suggested namespace
COBRAfields     = fieldFile{1}(dbFields);
fclose(fid);

% Load conversion between additional COBRA fields and namespaces:
fid             = fopen([ravenPath filesep 'struct_conversion' filesep 'cobraNamespaces.csv']);
fieldFile       = textscan(fid,'%s %s','Delimiter',',','HeaderLines',0);
COBRAfields     = [COBRAfields; fieldFile{1}];
COBRAnamespace  = [COBRAnamespace; fieldFile{2}];
rxnCOBRAfields  = COBRAfields(startsWith(COBRAfields,'rxn'));
rxnNamespaces   = COBRAnamespace(startsWith(COBRAfields,'rxn'));
metCOBRAfields  = COBRAfields(startsWith(COBRAfields,'met'));
metNamespaces   = COBRAnamespace(startsWith(COBRAfields,'met'));
geneCOBRAfields = COBRAfields(startsWith(COBRAfields,'gene'));
geneNamespaces  = COBRAnamespace(startsWith(COBRAfields,'gene'));
fclose(fid);

if isRaven
    %Firstly remove boundary metabolites
    model=simplifyModel(model);
end

% Keep fields that have identical names and content
newModel.S=model.S;
newModel.lb=model.lb;
newModel.ub=model.ub;
newModel.c=model.c;
newModel.rxns=model.rxns;
optFields = {'rxnNames','subSystems','rxnNotes',...
    'metFormulas','comps','compNames','metCharges','genes',...
    'rxnConfidenceScores','rxnGeneMat','metNotes','rev'};
for i=1:length(optFields)
    if isfield(model,optFields{i})
        newModel.(optFields{i})=model.(optFields{i});
    end
end
    
% Convert unique fields
if isRaven
    fprintf('Converting RAVEN structure to COBRA..\n');
    %Convert from RAVEN to COBRA structure
    
    %Mandatory COBRA fields
    newModel.rxns=model.rxns;
    if all(~cellfun(@isempty,regexp(model.mets,'\[[^\]]+\]$')))
        newModel.mets=model.mets;
    else
        %Check if model has compartment info as "met_c" suffix in all metabolites:
        BiGGformat = false(size(model.mets));
        for i=1:numel(model.comps)
            compPos=model.metComps==i;
            BiGGformat(compPos)=~cellfun(@isempty,regexp(model.mets(compPos),['_' model.comps{i} '$']));
        end
        if all(BiGGformat)
            newModel.mets=model.mets;
            for i=1:numel(model.comps)
                newModel.mets=regexprep(newModel.mets,['_' model.comps{i} '$'],['[' model.comps{i} ']']);
            end
        else
            newModel.mets=strcat(model.mets,'[',model.comps(model.metComps),']');
        end
    end

    %b, csense, osenseStr, genes, rules are also mandatory, but defined
    %later to match the order of fields
    
    %Optional COBRA fields
    if isfield(model,'id')
        newModel.modelID=model.id;
    end
    if isfield(model,'name')
        newModel.modelName=model.name;
    end
    if isfield(model,'eccodes')
        newModel.rxnECNumbers=model.eccodes;
    end
    if isfield(model,'rxnMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.rxnMiriams);
        for i = 1:length(rxnCOBRAfields)
            j=ismember(extractedMiriamNames,rxnNamespaces{i});
            if any(j)
                eval(['newModel.' rxnCOBRAfields{i} ' = miriams(:,j);'])
            end
        end
    end
    if isfield(model,'rxnReferences') % Concatenate model.rxnReferences to those extracted from model.rxnMiriams
        if isfield(newModel,'rxnReferences')
            newModel.rxnReferences = strcat(newModel.rxnReferences,{'; '},model.rxnReferences);
            newModel.rxnReferences = regexprep(newModel.rxnReferences,'^; $','');
        else
            newModel.rxnReferences = model.rxnReferences;
        end
    end
    if isfield(model,'metNames')
        newModel.metNames=strcat(model.metNames,' [',model.compNames(model.metComps),']');
    end
    if isfield(model,'metMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.metMiriams);
        %Shorten miriam names for KEGG and PubChem. These shorter names
        %will be used later to concatenate KEGG COMPOUND/GLYCAN and PubChem
        %Compound/Substance, into corresponding COBRA model fields
        extractedMiriamNames=regexprep(extractedMiriamNames,'^kegg\..+','kegg');
        extractedMiriamNames=regexprep(extractedMiriamNames,'^pubchem\..+','pubchem');
        i=ismember(extractedMiriamNames,'kegg');
        if any(i) % Combine KEGG compounds and glycans
            for j=1:length(i)
                if i(j) && isfield(newModel,'metKEGGID')~=1
                    newModel.metKEGGID=miriams(:,j);
                elseif i(j)
                    newModel.metKEGGID=strcat(newModel.metKEGGID,';',miriams(:,j));
                end
            end
            newModel.metKEGGID=regexprep(newModel.metKEGGID,'^;|;$','');
        end
        i=ismember(extractedMiriamNames,'pubchem');
        if any(i) % Combine Pubchem compounds and substances
            for j=1:length(i)
                if i(j) && isfield(newModel,'metPubChemID')~=1
                    newModel.metPubChemID=miriams(:,j);
                elseif i(j)
                    newModel.metPubChemID=strcat(newModel.metPubChemID,';',miriams(:,j));
                end
            end
            newModel.metPubChemID=regexprep(newModel.metPubChemID,'^;|;$','');
        end
        %All other Miriams can be directly parsed with no modifications:
        for i = 1:length(metCOBRAfields)
            j=ismember(extractedMiriamNames,metNamespaces{i});
            if any(j)
                eval(['newModel.' metCOBRAfields{i} ' = miriams(:,j);'])
            end
        end
    end
    if isfield(model,'inchis')
        newModel.metInChIString=regexprep(strcat('InChI=', model.inchis),'^InChI=$','');
    end
    newModel.b=zeros(numel(model.mets),1);
    newModel.csense=repmat('E',size(model.mets));
    if isfield(model,'geneMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.geneMiriams);
        for i = 1:length(geneCOBRAfields)
            j=ismember(extractedMiriamNames,geneNamespaces{i});
            if any(j)
                eval(['newModel.' geneCOBRAfields{i} ' = miriams(:,j);'])
            end
        end
    end
    if isfield(model,'geneShortNames')
        newModel.geneNames=model.geneShortNames;
    end
    if isfield(model,'genes')
        newModel.rules=grrulesToRules(model);
    else
        fprintf('WARNING: no genes detected. The model therefore may not be exportable to SBML file with writeCbModel\n');
    end
    newModel.osenseStr='max';
else
    fprintf('Converting COBRA structure to RAVEN..\n');
    %Convert from COBRA to RAVEN structure
    
    %Mandatory RAVEN fields
    newModel.mets=model.mets;
    if ~isfield(model,'comps')
        model.comps = unique(regexprep(model.mets,'.*\[([^\]]+)\]$','$1'));
    end
    for i=1:numel(model.comps)
        newModel.mets=regexprep(newModel.mets,['\[', model.comps{i}, '\]$'],'');
        newModel.mets=regexprep(newModel.mets,['\[', model.compNames{i}, '\]$'],'');
    end
    
    %In some cases (e.g. any model that uses BiGG ids as main ids), there
    %may be overlapping mets due to removal of compartment info. To avoid
    %this, we change compartments from e.g. [c] into _c
    if numel(unique(newModel.mets))~=numel(model.mets)
        newModel.mets=model.mets;
        for i=1:numel(model.comps)
            newModel.mets=regexprep(newModel.mets,['\[' model.comps{i} '\]$'],['_' model.comps{i}]);
        end
    end
    %Since COBRA no longer contains rev field it is assumed that rxn is
    %reversible if its lower bound is set below zero
    if ~isfield(model,'rev')
        for i=1:numel(model.rxns)
            if model.lb(i)<0
                newModel.rev(i,1)=1;
            else
                newModel.rev(i,1)=0;
            end
        end
    end
    newModel.b=zeros(numel(model.mets),1);
    if ~isfield(model,'comps')
        %Since 'comps' field is not mandatory in COBRA, it may be required
        %to obtain the non-redundant list of comps from metabolite ids, if
        %'comps' field is not available
        newModel.comps=regexprep(model.mets,'^.+\[','');
        newModel.comps=regexprep(newModel.comps,'\]$','');
        newModel.comps=unique(newModel.comps);
    end
    
    %metComps is also mandatory, but defined later to match the order of
    %fields
    
    %Fields 'name' and 'id' are also considered as mandatory, but
    %these are added to the model during exportModel/exportToExcelFormat
    %anyway, so there is no point to add this information here
    
    %Optional RAVEN fields
    if isfield(model,'modelID')
        newModel.id=model.modelID;
    end
    if isfield(model,'modelName')
        newModel.name=model.modelName;
    end
    if isfield(model,'rules')
        model.grRules        = rulesTogrrules(model);
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        newModel.grRules     = grRules;
        newModel.rxnGeneMat  = rxnGeneMat;
    end
    if isfield(model,'rxnECNumbers')
        newModel.eccodes=regexprep(model.rxnECNumbers,'EC|EC:','');
    end
    if any(isfield(model,rxnCOBRAfields))
        for i=1:numel(model.rxns)
            counter=1;
            newModel.rxnMiriams{i,1}=[];
            if isfield(model,'rxnReferences')
                if ~isempty(model.rxnReferences{i})
                    pmids = model.rxnReferences{i};
                    pmids = strsplit(pmids,'; ');
                    nonPmids = cellfun(@isempty,regexp(pmids,'^\d+$','match','once'));
                    if any(nonPmids) %Not a pubmed id, keep in rxnReferences instead
                        newModel.rxnReferences{i,1} = strjoin(pmids(nonPmids),', ');
                        pmids(nonPmids)=[];
                    end
                    for j = 1:length(pmids)
                        newModel.rxnMiriams{i,1}.name{counter,1} = 'pubmed';
                        newModel.rxnMiriams{i,1}.value{counter,1} = pmids{j};
                        counter=counter+1;
                    end
                end
            end
            for j = 2:length(rxnCOBRAfields) %Start from 2, as 1 is rxnReferences
                if isfield(model,rxnCOBRAfields{j})
                    rxnAnnotation = eval(['model.' rxnCOBRAfields{j} '{i}']);
                    if ~isempty(rxnAnnotation)
                        rxnAnnotation = strtrim(strsplit(rxnAnnotation,';'));
                        for a=1:length(rxnAnnotation)
                            newModel.rxnMiriams{i,1}.name{counter,1} = rxnNamespaces{j};
                            newModel.rxnMiriams{i,1}.value{counter,1} = rxnAnnotation{a};
                            counter=counter+1;
                        end
                    end
                end
            end
        end
    end
    if isfield(newModel,'rxnReferences')
        emptyEntry = cellfun(@isempty,newModel.rxnReferences);
        newModel.rxnReferences(emptyEntry)={''};
    end
    if any(isfield(model,geneCOBRAfields))
        for i=1:numel(model.genes)
            counter=1;
            newModel.geneMiriams{i,1}=[];
            for j = 1:length(geneCOBRAfields)
                if isfield(model,geneCOBRAfields{j})
                    geneAnnotation = eval(['model.' geneCOBRAfields{j} '{i}']);
                    if ~isempty(geneAnnotation)
                        geneAnnotation = strtrim(strsplit(geneAnnotation,';'));
                        for a=1:length(geneAnnotation)
                            newModel.geneMiriams{i,1}.name{counter,1} = geneNamespaces{j};
                            newModel.geneMiriams{i,1}.value{counter,1} = geneAnnotation{a};
                            counter=counter+1;
                        end
                    end
                end
            end
        end
    end
    if isfield(model,'geneNames')
        newModel.geneShortNames=model.geneNames;
    end
    newModel.metNames=model.metNames;
    for i=1:numel(model.comps)
        newModel.metNames=regexprep(newModel.metNames,['\[', model.comps{i}, '\]$'],'');
        newModel.metNames=regexprep(newModel.metNames,['\[', model.compNames{i}, '\]$'],'');
    end
    newModel.metNames=deblank(newModel.metNames);
    newModel.metComps=regexprep(model.mets,'^.+\[','');
    newModel.metComps=regexprep(newModel.metComps,'\]$','');
    [~, newModel.metComps]=ismember(newModel.metComps,newModel.comps);
    if isfield(model,'metInChIString')
        newModel.inchis=regexprep(model.metInChIString,'^InChI=','');
    end
    printWarning=false;
    if any(isfield(model,[metCOBRAfields;'metKEGGID';'metPubChemID']))
        for i=1:numel(model.mets)
            counter=1;
            newModel.metMiriams{i,1}=[];
            if isfield(model,'metKEGGID')
                if ~isempty(model.metKEGGID{i})
                    if strcmp(model.metKEGGID{i}(1),'C')
                        newModel.metMiriams{i,1}.name{counter,1} = 'kegg.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    elseif strcmp(model.metKEGGID{i}(1),'G')
                        newModel.metMiriams{i,1}.name{counter,1} = 'kegg.glycan';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    end
                end
            end
            if isfield(model,'metPubChemID')
                if ~isempty(model.metPubChemID{i})
                    if length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'CID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    elseif length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'SID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.substance';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    else
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                        printWarning=true;
                    end
                end
            end            
            for j = 1:length(metCOBRAfields)
                if isfield(model,metCOBRAfields{j})
                    metAnnotation = eval(['model.' metCOBRAfields{j} '{i}']);
                    if ~isempty(metAnnotation)
                        metAnnotation = strtrim(strsplit(metAnnotation,';'));
                        for a=1:length(metAnnotation)
                            newModel.metMiriams{i,1}.name{counter,1} = metNamespaces{j};
                            newModel.metMiriams{i,1}.value{counter,1} = metAnnotation{a};
                            counter=counter+1;
                        end
                    end
                end
            end
        end
    end
    if printWarning
        fprintf('Could not determine whether PubChemIDs are compounds (CID)\n or substances (SID). All annotated PubChemIDs will therefore \n be assigned as compounds (CID).\n');
    end
end

% Order fields
modelNew=standardizeModelFieldOrder(newModel); % Corrects for both RAVEN and COBRA models
end

function rules=grrulesToRules(model)
%This function just takes grRules, changes all gene names to
%'x(geneNumber)' and also changes 'or' and 'and' relations to corresponding
%symbols
replacingGenes=cell([size(model.genes,1) 1]);
for i=1:numel(replacingGenes)
    replacingGenes{i}=strcat('x(',num2str(i),')');
end
rules = strcat({' '},model.grRules,{' '});
for i=1:length(model.genes)
    rules=regexprep(rules,[' ' model.genes{i} ' '],[' ' replacingGenes{i} ' ']);
    rules=regexprep(rules,['(' model.genes{i} ' '],['(' replacingGenes{i} ' ']);
    rules=regexprep(rules,[' ' model.genes{i} ')'],[' ' replacingGenes{i} ')']);
end
rules=regexprep(rules,' and ',' & ');
rules=regexprep(rules,' or ',' | ');
rules=strtrim(rules);
end

function grRules=rulesTogrrules(model)
%This function takes rules, replaces &/| for and/or, replaces the x(i)
%format with the actual gene ID, and takes out extra whitespace and
%redundant parenthesis introduced by COBRA, to create grRules.
grRules = strrep(model.rules,'&','and');
grRules = strrep(grRules,'|','or');
for i = 1:length(model.genes)
    grRules = strrep(grRules,['x(' num2str(i) ')'],model.genes{i});
end
grRules = strrep(grRules,'( ','(');
grRules = strrep(grRules,' )',')');
grRules = regexprep(grRules,'^(',''); %rules that start with a "("
grRules = regexprep(grRules,')$',''); %rules that end with a ")"
end
