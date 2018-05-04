function model=importExcelModel(fileName,removeExcMets,printWarnings,ignoreErrors)
% importExcelModel
%   Imports a constraint-based model from a Excel file
%
%   fileName      a Microsoft Excel file to import
%   removeExcMets true if exchange metabolites should be removed. This is
%                 needed to be able to run simulations, but it could also
%                 be done using simplifyModel at a later stage (opt,
%                 default true)
%   printWarnings true if warnings should be printed (opt, default true)
%   ignoreErrors  true if errors should be ignored. See below for details
%                 (opt, default false)
%
%   model
%       description      description of model contents
%       id               model ID
%       rxns             reaction ids
%       mets             metabolite ids
%       S                stoichiometric matrix
%       lb               lower bounds
%       ub               upper bounds
%       rev              reversibility vector
%       c                objective coefficients
%       b                equality constraints for the metabolite equations
%       comps            compartment ids
%       compNames        compartment names
%       compOutside      the id (as in comps) for the compartment
%                        surrounding each of the compartments
%       compMiriams      structure with MIRIAM information about the
%                        compartments
%       rxnNames         reaction description
%       rxnComps         compartments for reactions
%       grRules          reaction to gene rules in text form
%       rxnGeneMat       reaction-to-gene mapping in sparse matrix form
%       subSystems       subsystem name for each reaction
%       eccodes          EC-codes for the reactions
%       rxnMiriams       structure with MIRIAM information about the reactions
%       rxnNotes         reaction notes
%       rxnReferences    reaction references
%       rxnConfidenceScores reaction confidence scores
%       genes            list of all genes
%       geneComps        compartments for genes
%       geneMiriams      structure with MIRIAM information about the genes
%       geneShortNames   gene alternative names (e.g. ERG10)
%       metNames         metabolite description
%       metComps         compartments for metabolites
%       inchis           InChI-codes for metabolites
%       metFormulas      metabolite chemical formula
%       metMiriams       structure with MIRIAM information about the metabolites
%       metCharges       metabolite charge
%       unconstrained    true if the metabolite is an exchange metabolite
%
%   Loads models in the RAVEN Toolbox Excel format. A number of consistency
%   checks are performed in order to ensure that the model is valid. These
%   can be ignored by putting ignoreErrors to true. However, this is highly
%   advised against, as it can result in errors in simulations or other
%   functionalities. The RAVEN Toolbox is made to function only on consistent
%   models, and the only checks performed are when the model is imported.
%
%   NOTE: Most errors are checked for by checkModelStruct, but some
%   are checked for in this function as well. Those are ones which relate
%   to missing model elements and so on, and which would make it impossible
%   to construct the model structure. Those errors cannot be ignored by
%   setting ignoreErrors to true.
%
%   Usage: model=importExcelModel(fileName,removeExcMets,printWarnings,ignoreErrors)
%
%   Simonas Marcisauskas, 2018-04-04
%

if nargin<2
    removeExcMets=true;
end
if nargin<3
    printWarnings=true;
end
if nargin<4
    ignoreErrors=false;
end

%This is to match the order of the fields to those you get from importing
%from SBML
model=[];
model.id=[];
model.description=[];
model.annotation=[];
%Default bounds if not defined
model.annotation.defaultLB=-1000;
model.annotation.defaultUB=1000;
model.rxns={};
model.mets={};
model.S=[];
model.lb=[];
model.ub=[];
model.rev=[];
model.c=[];
model.b=[];
model.comps={};
model.compNames={};
model.compOutside={};
model.compMiriams={};
model.rxnNames={};
model.rxnComps={}; %Will be double later
model.grRules={};
model.rxnGeneMat=[];
model.subSystems={};
model.eccodes={};
model.rxnMiriams={};
model.rxnNotes={};
model.rxnReferences={};
model.rxnConfidenceScores={}; %Will be double later
model.genes={};
model.geneComps={}; %Will be double later
model.geneMiriams={};
model.geneShortNames={};
model.metNames={};
model.metComps=[];
model.inchis={};
model.metFormulas={};
model.metMiriams={};
model.metCharges={}; %Will be double later
model.unconstrained=[];

workbook=loadWorkbook(fileName);

[raw, flag]=loadSheet(workbook,'MODEL');

if flag<0
    if printWarnings==true
        EM='Could not load the MODEL sheet';
        dispEM(EM,false);
    end
    model.id='UNKNOWN';
    model.description='No model details available';
else
    raw=cleanSheet(raw);
    
    %It is assumed that the first line is labels and that the second one is
    %info
    allLabels={'ID';'DESCRIPTION';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'TAXONOMY';'NOTES'};
    
    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'MODELID','ID');
    raw(1,:)=strrep(raw(1,:),'MODELNAME','DESCRIPTION');
    
    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
                if any(raw{I(i),2})
                    model.id=toStr(raw{2,I(i)}); %Should be string already
                else
                    EM='No model ID supplied';
                    dispEM(EM);
                end
            case 2
                if any(raw{2,I(i)})
                    model.description=toStr(raw{2,I(i)}); %Should be string already
                else
                    EM='No model name supplied';
                    dispEM(EM);
                end
            case 3
                if ~isempty(raw{2,I(i)})
                    try
                        model.annotation.defaultLB=toDouble(raw{2,I(i)},NaN);
                    catch
                        EM='DEFAULT LOWER must be numeric';
                        dispEM(EM);
                    end
                else
                    if printWarnings==true
                        fprintf('NOTE: DEFAULT LOWER not supplied. Uses -1000\n');
                    end
                    model.annotation.defaultLB=-1000;
                end
            case 4
                if ~isempty(raw{2,I(i)})
                    try
                        model.annotation.defaultUB=toDouble(raw{2,I(i)},NaN);
                    catch
                        EM='DEFAULT UPPER must be numeric';
                        dispEM(EM);
                    end
                else
                    if printWarnings==true
                        fprintf('NOTE: DEFAULT UPPER not supplied. Uses 1000\n');
                    end
                    model.annotation.defaultUB=1000;
                end
            case 5
                if any(raw{2,I(i)})
                    model.annotation.givenName=toStr(raw{2,I(i)}); %Should be string already
                end
            case 6
                if any(raw{2,I(i)})
                    model.annotation.familyName=toStr(raw{2,I(i)}); %Should be string already
                end
            case 7
                if any(raw{2,I(i)})
                    model.annotation.email=toStr(raw{2,I(i)}); %Should be string already
                end
            case 8
                if any(raw{2,I(i)})
                    model.annotation.organization=toStr(raw{2,I(i)}); %Should be string already
                end
            case 9
                if any(raw{2,I(i)})
                    model.annotation.taxonomy=toStr(raw{2,I(i)}); %Should be string already
                end
            case 10
                if any(raw{2,I(i)})
                    model.annotation.note=toStr(raw{2,I(i)}); %Should be string already
                end
        end
    end
end

%Get compartment information
[raw, flag]=loadSheet(workbook,'COMPS');

if flag<0
    if printWarnings==true
        EM='Could not load the COMPS sheet. All elements will be assigned to a compartment "s" for "System"';
        dispEM(EM,false);
    end
    model.comps={'s'};
    model.compNames={'System'};
else
    raw=cleanSheet(raw);
    
    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'COMPABBREV','ABBREVIATION');
    raw(1,:)=strrep(raw(1,:),'COMPNAME','NAME');
    
    allLabels={'ABBREVIATION';'NAME';'INSIDE';'MIRIAM'};
    
    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
                model.comps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 2
                model.compNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                model.compOutside=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 4
                model.compMiriams=raw(2:end,I(i));
        end
    end
    
    %Check that necessary fields are loaded
    if isempty(model.comps)
        EM='There must be a column named ABBREVIATION in the COMPS sheet';
        dispEM(EM);
    end
    if isempty(model.compNames)
        model.compNames=model.comps;
        if printWarnings==true
            EM='There is no column named NAME in the COMPS sheet. ABBREVIATION will be used as name';
            dispEM(EM,false);
        end
    end
    
    model.compMiriams=parseMiriam(model.compMiriams);
end

%Get all the genes and info about them
[raw, flag]=loadSheet(workbook,'GENES');

if flag<0
    if printWarnings==true
        EM='There is no spreadsheet named GENES';
        dispEM(EM,false)
    end
else
    raw=cleanSheet(raw);
    
    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'GENE NAME','NAME');
    
    allLabels={'NAME';'MIRIAM';'SHORT NAME';'COMPARTMENT'};
    
    %Loop through the labels
    [I, J]=ismember(upper(raw(1,:)),allLabels);
    I=find(I);
    foundGenes=false;
    for i=1:numel(I)
        switch J(I(i))
            case 1
                model.genes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
                foundGenes=true;
            case 2
                model.geneMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                model.geneShortNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 4
                model.geneComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        end
    end
    
    if foundGenes==false
        EM='There must be a column named NAME in the GENES sheet';
        dispEM(EM);
    end
    
    %Its ok if all of them are empty
    if all(cellfun(@isempty,model.geneComps))
        model.geneComps=[];
    end
    
    %Check that geneName contain only strings and no empty strings
    if ~iscellstr(model.genes)
        EM='All gene names have to be strings';
        dispEM(EM);
    else
        if any(strcmp('',model.genes))
            EM='There can be no empty strings in gene names';
            dispEM(EM);
        end
    end
    
    %Check that geneComp contains only strings and no empty string
    if ~isempty(model.geneComps)
        if ~iscellstr(model.geneComps)
            EM='All gene compartments have to be strings';
            dispEM(EM);
        else
            if any(strcmp('',model.geneComps))
                EM='There can be no empty strings in gene compartments';
                dispEM(EM);
            end
        end
        [I, model.geneComps]=ismember(model.geneComps,model.comps);
        EM='The following genes have compartment abbreviations which could not be found:';
        dispEM(EM,true,model.genes(~I));
    end
end

model.geneMiriams=parseMiriam(model.geneMiriams);

%Loads the reaction data
[raw, flag]=loadSheet(workbook,'RXNS');

if flag<0
    EM='Could not load the RXNS sheet';
    dispEM(EM);
end

raw=cleanSheet(raw);

%Map to new captions
raw(1,:)=upper(raw(1,:));
raw(1,:)=strrep(raw(1,:),'RXNID','ID');

allLabels={'ID';'NAME';'EQUATION';'EC-NUMBER';'GENE ASSOCIATION';'LOWER BOUND';'UPPER BOUND';'OBJECTIVE';'COMPARTMENT';'SUBSYSTEM';'REPLACEMENT ID';'MIRIAM';'NOTE';'REFERENCE';'CONFIDENCE SCORE'};
equations={};
reactionReplacement={};

%Loop through the labels
[I, J]=ismember(raw(1,:),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
            model.rxns=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 2
            model.rxnNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 3
            equations=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 4
            model.eccodes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 5
            model.grRules=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 6
            try
                model.lb=cellfun(@(x) toDouble(x,NaN),raw(2:end,I(i)));
            catch
                EM='The lower bounds must be numerical values';
                dispEM(EM);
            end
        case 7
            try
                model.ub=cellfun(@(x) toDouble(x,NaN),raw(2:end,I(i)));
            catch
                EM='The upper bounds must be numerical values';
                dispEM(EM);
            end
        case 8
            try
                model.c=cellfun(@(x) toDouble(x,0),raw(2:end,I(i)));
            catch
                EM='The objective coefficients must be numerical values';
                dispEM(EM);
            end
        case 9
            model.rxnComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 10
            subsystems=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 11
            reactionReplacement=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 12
            model.rxnMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 13
            model.rxnNotes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 14
            model.rxnReferences=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 15
            model.rxnConfidenceScores=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
    end
end

if ~isempty(model.rxnConfidenceScores)
    model.rxnConfidenceScores=str2double(model.rxnConfidenceScores);
end
for i=1:numel(subsystems)
    model.subSystems{i,1}=cellstr(strsplit(subsystems{i,1},';'));
end

%Check that all necessary reaction info has been loaded
if isempty(equations)
    EM='There must be a column named EQUATION in the RXNS sheet';
    dispEM(EM);
end
if isempty(model.rxns)
    if printWarnings==true
        EM='There is no column named ID in the RXNS sheet. The reactions will be named as "r1", "r2"...';
        dispEM(EM,false);
    end
    I=num2cell((1:numel(equations))');
    model.rxns=strcat('r',cellfun(@num2str,I,'UniformOutput',false));
end

%Check if some other stuff is loaded and populate with default values
%otherwise
if isempty(model.rxnNames)
    model.rxnNames=cell(numel(model.rxns),1);
    model.rxnNames(:)={''};
    if printWarnings==true
        EM='There is no column named NAME in the RXNS sheet. Empty strings will be used as reaction names';
        dispEM(EM,false);
    end
end
if isempty(model.lb)
    %This is not set here since the reversibility isn't known yet
    model.lb=nan(numel(model.rxns),1);
    if printWarnings==true
        EM='There is no column named LOWER BOUND in the RXNS sheet. Default bounds will be used';
        dispEM(EM,false);
    end
end
if isempty(model.ub)
    model.ub=nan(numel(model.rxns),1);
    if printWarnings==true
        EM='There is no column named UPPER BOUND in the RXNS sheet. Default bounds will be used';
        dispEM(EM,false);
    end
end
if isempty(model.c)
    model.c=zeros(numel(model.rxns),1);
    if printWarnings==true
        EM='There is no column named OBJECTIVE in the RXNS sheet';
        dispEM(EM,false);
    end
end

%Either all reactions must have a compartment string or none of them. Check
%if it's only empty and if so return it to []
if ~isempty(model.rxnComps)
    if all(cellfun(@isempty,model.rxnComps))
        model.rxnComps=[];
    end
end

%Construct the rxnMiriams structure
model.rxnMiriams=parseMiriam(model.rxnMiriams);

%Replace the reaction IDs for those IDs that have a corresponding
%replacement name.
I=cellfun(@any,reactionReplacement);
model.rxns(I)=reactionReplacement(I);

%Check that there are no empty strings in reactionIDs or equations
if any(strcmp('',model.rxns))
    EM='There are empty reaction IDs';
    dispEM(EM);
end

if any(strcmp('',equations))
    EM='There are empty equations';
    dispEM(EM);
end

if ~isempty(model.rxnComps)
    if any(strcmp('',model.rxnComps))
        EM='Either all reactions must have an associated compartment string or none of them';
        dispEM(EM);
    end
end

if ~isempty(model.grRules)
    tempRules=model.grRules;
    for i=1:length(model.rxns)
        %Check that all gene associations have a match in the gene list
        if ~isempty(model.grRules{i})
            tempRules{i}=regexprep(tempRules{i},' and | or ','>'); %New format: Genes are separated 'and' and 'or' strings with parentheses
            tempRules{i}=regexprep(tempRules{i},'(',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
            tempRules{i}=regexprep(tempRules{i},')',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
            indexesNew=strfind(tempRules{i},'>'); %Old format: Genes are separated by ":" for AND and ";" for OR
            indexes=strfind(tempRules{i},':'); %Old format: Genes are separated by ":" for AND and ";" for OR
            indexes=unique([indexesNew indexes strfind(tempRules{i},';')]);
            if isempty(indexes)
                %See if you have a match
                I=find(strcmp(tempRules{i},model.genes));
                if isempty(I)
                    EM=['The gene association in reaction ' model.rxns{i} ' (' tempRules{i} ') is not present in the gene list'];
                    dispEM(EM);
                end
            else
                temp=[0 indexes numel(tempRules{i})+1];
                for j=1:numel(indexes)+1
                    %The reaction has several associated genes
                    geneName=tempRules{i}(temp(j)+1:temp(j+1)-1);
                    I=find(strcmp(geneName,model.genes));
                    if isempty(I)
                        EM=['The gene association in reaction ' model.rxns{i} ' (' geneName ') is not present in the gene list'];
                        dispEM(EM);
                    end
                end
            end
            %In order to adhere to the COBRA standards it should be like
            %this: -If only one gene then no parentheses -If only "and" or
            %only "or" there should only be one set of parentheses -If both
            %"and" and "or", then split on "or". This is not complete, but
            %it's the type of relationship supported by the Excel
            %formulation
            aSign=strfind(model.grRules{i},':');
            oSign=strfind(model.grRules{i},';');
            if isempty(aSign) && isempty(oSign)
                model.grRules{i}=model.grRules{i};
            else
                if isempty(aSign)
                    model.grRules{i}=['(' strrep(model.grRules{i},';',' or ') ')'];
                else
                    if isempty(oSign)
                        model.grRules{i}=['(' strrep(model.grRules{i},':',' and ') ')'];
                    else
                        model.grRules{i}=['((' strrep(model.grRules{i},';',') or (') '))'];
                        model.grRules{i}=strrep(model.grRules{i},':',' and ');
                    end
                end
            end
        end
    end
end

%Check that the compartment for each reaction can be found
if ~isempty(model.rxnComps)
    [I, model.rxnComps]=ismember(model.rxnComps,model.comps);
    EM='The following reactions have compartment abbreviations which could not be found:';
    dispEM(EM,true,model.rxns(~I));
end

%Get all the metabolites and info about them
[raw, flag]=loadSheet(workbook,'METS');

if flag<0
    if printWarnings==true
        EM='There is no spreadsheet named METS. The metabolites will be named "m1", "m2"... and assigned to the first compartment';
        dispEM(EM,false);
    end
    %Parse the equations to find out how many metabolites there are
    metsForParsing=parseRxnEqu(equations);
    I=num2cell((1:numel(metsForParsing))');
    model.mets=strcat('m',cellfun(@num2str,I,'UniformOutput',false));
    model.metComps=ones(numel(model.mets),1);
    model.unconstrained=zeros(numel(model.mets),1);
    model.metNames=metsForParsing;
else
    raw=cleanSheet(raw);
    
    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'METID','ID');
    raw(1,:)=strrep(raw(1,:),'METNAME','NAME');
    
    allLabels={'ID';'NAME';'UNCONSTRAINED';'MIRIAM';'COMPOSITION';'INCHI';'COMPARTMENT';'REPLACEMENT ID';'CHARGE'};
    
    %Load the metabolite information
    metReplacement={};
    
    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
                model.mets=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 2
                model.metNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                model.unconstrained=cellfun(@boolToDouble,raw(2:end,I(i)));
                
                %NaN is returned if the values couldn't be parsed
                EM='The UNCONSTRAINED property for the following metabolites must be "true"/"false", 1/0, TRUE/FALSE or not set:';
                dispEM(EM,true,model.mets(isnan(model.unconstrained)));
            case 4
                model.metMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 5
                model.metFormulas=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 6
                model.inchis=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 7
                model.metComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
                
                %Check that all metabolites have compartments defined
                if any(strcmp('',model.metComps))
                    EM='All metabolites must have an associated compartment string';
                    dispEM(EM);
                end
            case 8
                metReplacement=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 9
                model.metCharges=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        end
    end
    
    %Check that necessary fields are loaded (METID)
    if isempty(model.mets)
        EM='There must be a column named ID in the METS sheet';
        dispEM(EM);
    end
    
    %Check that some other stuff is loaded and use default values otherwise
    if isempty(model.metNames)
        model.metNames=cell(numel(model.mets),1);
        if printWarnings==true
            EM='There is no column named NAME in the METS sheet. ID will be used as name';
            dispEM(EM,false);
        end
    end
    if isempty(model.unconstrained)
        model.unconstrained=zeros(numel(model.mets),1);
        if printWarnings==true
            EM='There is no column named UNCONSTRAINED in the METS sheet. All metabolites will be constrained';
            dispEM(EM,false);
        end
    end
    
    if isempty(model.metComps)
        model.metComps=cell(numel(model.mets),1);
        model.metComps(:)=model.comps(1);
        if printWarnings==true
            EM='There is no column named COMPARTMENT in the METS sheet. All metabolites will be assigned to the first compartment in COMPS. Note that RAVEN makes extensive use of metabolite names and compartments. Some features will therefore not function correctly if metabolite compartments are not correctly assigned';
            dispEM(EM,false);
        end
    end
    
    %The composition should be loaded from InChIs when available
    I=find(~cellfun(@isempty,model.inchis));
    for i=1:numel(I)
        S=regexp(model.inchis(I(i)),'/','split');
        S=S{1};
        if numel(S)>=2
            %Don't copy if it doesn't look good
            model.metFormulas(I(i))=S(2);
        end
    end
    
    %Check that the compartment for each metabolite can be found. Also
    %convert from id to index
    [I, model.metComps]=ismember(model.metComps,model.comps);
    EM='The following metabolites have compartment abbreviations which could not be found:';
    dispEM(EM,true,model.mets(~I));
    
    %Check that the model.mets vector is unique. The problem is that the
    %checkModelStruct cannot check for that since only the metReplacements
    %(if used) end up in the model structure
    I=false(numel(model.mets),1);
    [J, K]=unique(model.mets);
    if numel(J)~=numel(model.mets)
        L=1:numel(model.mets);
        L(K)=[];
        I(L)=true;
    end
    EM='The following metabolites are duplicates:';
    dispEM(EM,~ignoreErrors,model.mets(I));
    
    %Check that there are no metabolite IDs which are numbers. This would
    %give errors when parsing the equations
    I=cellfun(@str2double,model.mets);
    EM='The following metabolites have names which cannot be distinguished from numbers:';
    dispEM(EM,~ignoreErrors,model.mets(~isnan(I)));
    I=cellfun(@str2double,metReplacement);
    EM='The following metabolites have names which cannot be distinguished from numbers:';
    dispEM(EM,~ignoreErrors,metReplacement(~isnan(I)));
    
    %Replace the metabolite IDs for those IDs that have a corresponding
    %replacement metabolite. This is not used for matching, but will be
    %checked for consistency with SBML naming conventions
    metsForParsing=model.mets; %This is because the equations are written with this
    I=cellfun(@any,metReplacement);
    model.mets(I)=metReplacement(I);
    
    %If the metabolite name isn't set, replace it with the metabolite id
    I=~cellfun(@any,model.metNames);
    model.metNames(I)=model.mets(I);
    
    %Construct the metMiriams structure
    model.metMiriams=parseMiriam(model.metMiriams);
    
    %Either all metabolites have charge or none of them. Check if it's only
    %empty and if so return it to []
    if ~isempty(model.metCharges)
        if all(cellfun(@isempty,model.metCharges))
            model.metCharges=[];
        end
    end
    if ~isempty(model.metCharges)
        model.metCharges=str2double(model.metCharges);
    end
end

%Everything seems fine with the metabolite IDs, compartments, genes, and
%reactions

%Parse the equations
[model.S, mets, badRxns, model.rev]=constructS(equations,metsForParsing,model.rxns);
model.rev=model.rev*1; %Typecast to double

%Add default constraints
model.lb(isnan(model.lb))=model.annotation.defaultLB.*model.rev(isnan(model.lb));
model.ub(isnan(model.ub))=model.annotation.defaultUB;

%Reorder the S matrix so that it fits with the metabolite list in the
%structure
[~, I]=ismember(mets,metsForParsing);
model.S=model.S(I,:);

%Print warnings about the reactions which contain the same metabolite as
%both reactants and products
EM='The following reactions have metabolites which are present more than once. Only the net reactions will be exported:';
dispEM(EM,false,model.rxns(badRxns));

model.b=zeros(numel(model.mets),1);

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;

%Remove unused fields
if all(cellfun(@isempty,model.compOutside))
    model=rmfield(model,'compOutside');
end
if all(cellfun(@isempty,model.compMiriams))
    model=rmfield(model,'compMiriams');
end
if all(cellfun(@isempty,model.rxnNames))
    model=rmfield(model,'rxnNames');
end
if isempty(model.rxnComps)
    model=rmfield(model,'rxnComps');
end
if all(cellfun(@isempty,model.grRules))
    model=rmfield(model,'grRules');
end
if isfield(model,'rxnGeneMat') && isempty(model.rxnGeneMat)
    model=rmfield(model,'rxnGeneMat');
end
if all(cellfun(@isempty,model.subSystems))
    model=rmfield(model,'subSystems');
end
if all(cellfun(@isempty,model.eccodes))
    model=rmfield(model,'eccodes');
end
if all(cellfun(@isempty,model.rxnMiriams))
    model=rmfield(model,'rxnMiriams');
end
if all(cellfun(@isempty,model.rxnNotes))
    model=rmfield(model,'rxnNotes');
end
if all(cellfun(@isempty,model.rxnReferences))
    model=rmfield(model,'rxnReferences');
end
if isempty(model.rxnConfidenceScores)
    model=rmfield(model,'rxnConfidenceScores');
end
if isempty(model.genes)
    model=rmfield(model,'genes');
end
if isempty(model.geneComps)
    model=rmfield(model,'geneComps');
end
if isempty(model.geneMiriams)
    model=rmfield(model,'geneMiriams');
end
if all(cellfun(@isempty,model.geneShortNames))
    model=rmfield(model,'geneShortNames');
end
if all(cellfun(@isempty,model.inchis))
    model=rmfield(model,'inchis');
end
if all(cellfun(@isempty,model.metFormulas))
    model=rmfield(model,'metFormulas');
end
if all(cellfun(@isempty,model.metMiriams))
    model=rmfield(model,'metMiriams');
end
if isempty(model.metCharges)
    model=rmfield(model,'metCharges');
end

%The model structure has now been reconstructed but it can still contain
%many types of errors. The checkModelConsistency function is used to make
%sure that naming and mapping of stuff looks good
checkModelStruct(model,~ignoreErrors);

if removeExcMets==true
    model=simplifyModel(model);
end
end

function miriamStruct=parseMiriam(strings,miriamStruct)
%Gets the names and values of Miriam-string. Nothing fancy at all, just to
%prevent using the same code for metabolites, genes, and reactions. The
%function also allows for supplying a miriamStruct and the info will then
%be added

if nargin<2
    miriamStruct=cell(numel(strings),1);
end
for i=1:numel(strings)
    if any(strings{i})
        %A Miriam string can be several ids separated by ";". Each id is
        %"name(..:..)/value"; an old format when value is separated by
        %colon is also supported
        I=regexp(strings{i},';','split');
        if isfield(miriamStruct{i},'name')
            startIndex=numel(miriamStruct{i}.name);
            miriamStruct{i}.name=[miriamStruct{i}.name;cell(numel(I),1)];
            miriamStruct{i}.value=[miriamStruct{i}.value;cell(numel(I),1)];
        else
            startIndex=0;
            miriamStruct{i}.name=cell(numel(I),1);
            miriamStruct{i}.value=cell(numel(I),1);
        end
        
        for j=1:numel(I)
            if any(strfind(I{j},'/'))
                index=max(strfind(I{j},'/'));
            elseif any(strfind(I{j},':'))
                index=max(strfind(I{j},':'));
            end
            if any(index)
                miriamStruct{i}.name{startIndex+j}=I{j}(1:index-1);
                miriamStruct{i}.value{startIndex+j}=I{j}(index+1:end);
            else
                EM=['"' I{j} '" is not a valid MIRIAM string. The format must be "identifier/value" or identifier:value'];
                dispEM(EM);
            end
        end
    end
end
end

%For converting a value to string. This is used instead of num2str because
%I want to convert empty cells to {''}.
function y=toStr(x)
%x can be empty, numerical, string or boolean. It cannot be NaN. Boolean
%values will be converted to '1'/'0'
if isempty(x)
    y='';
else
    y=num2str(x);
end
end

%For converting to numeric. This is used instead of str2num because I want
%to be able to choose what empty values should be mapped to.
%
% default the value to use for empty input
function y=toDouble(x,default)
if isempty(x) %Note that this catches '' as well
    y=default;
else
    if isnumeric(x)
        y=x;
    else
        y=str2double(x);
        
        %This happens if the input couldn't be converted. Note that the
        %input itself cannot be NaN since it was fixed in clean imported
        if isnan(y)
            EM=['Cannot convert the string "' x '" to double'];
            dispEM(EM);
        end
    end
end
end

%For converting boolean (the UNCONSTRAINED field) to double (the
%model.unconstrained field)
function y=boolToDouble(x)
if isempty(x)
    y=0;
    return;
end
if islogical(x)
    y=x*1; %Typecast to double
    return;
end
if isnumeric(x)
    if x~=0
        y=1;
        return;
    else
        y=0;
        return;
    end
end
if ischar(x)
    if strcmpi(x,'TRUE')
        y=1;
        return;
    end
    if strcmpi(x,'FALSE')
        y=0;
        return;
    end
end
y=NaN; %This means that the input couldn't be parsed
end
