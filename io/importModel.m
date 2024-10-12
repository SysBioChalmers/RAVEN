function model=importModel(fileName,removeExcMets,COBRAstyle,supressWarnings)
% importModel
%   Import a constraint-based model from a SBML file
%
% Input:
%   fileName        a SBML file to import. A dialog window will open if
%                   no file name is specified.
%   removeExcMets   true if exchange metabolites should be removed. This is
%                   needed to be able to run simulations, but it could also
%                   be done using simplifyModel at a later stage (optional,
%                   default true)
%   COBRAstyle      true if the model uses COBRA style identifier prefixes
%                   that should be removed when loading the model: G_ for
%                   genes, R_ for reactions, M_ for metabolites, and C_ for
%                   compartments. (optional, default false)
%   supressWarnings true if warnings regarding the model structure should
%                   be supressed (optional, default false)
%
% Output:
%   model
%       id               model ID
%       name             name of model contents
%       annotation       additional information about model
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
%       proteins     protein associated to each gene
%       metNames         metabolite description
%       metComps         compartments for metabolites
%       inchis           InChI-codes for metabolites
%       metFormulas      metabolite chemical formula
%       metMiriams       structure with MIRIAM information about the metabolites
%       metCharges       metabolite charge
%       unconstrained    true if the metabolite is an exchange metabolite
%
% A number of consistency checks are performed in order to ensure that the
% model is valid. Take these warnings seriously and modify the model
% structure to solve them.
%
% Usage: model = importModel(fileName, removeExcMets, COBRAstyle, supressWarnings)

if nargin<1 || isempty(fileName)
    [fileName, pathName] = uigetfile({'*.xml;*.sbml'}, 'Please select the model file');
    if fileName == 0
        error('You should select a model file')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);
if nargin<2 || isempty(removeExcMets)
    removeExcMets=true;
end

if nargin<3 || isempty(COBRAstyle)
    COBRAstyle=false;
end

if nargin<4
    supressWarnings=false;
end

fileName=checkFileExistence(fileName,1);
% If path contains non-ASCII characters, copy file to tempdir first, as
% libSBML is known to have problems with this on Windows:
% https://sbml.org/software/libsbml/libsbml-docs/known-pitfalls/#matlab-on-windows-has-issues-with-unicode-filenames
if ispc && any(double(fileName)>128)
    [~,originalFile,ext] = fileparts(fileName);
    tempFile = fullfile(tempdir,[originalFile ext]);
    copyfile(fileName,tempFile);
    fileName = tempFile;
end

%This is to match the order of the fields to those you get from importing
%from Excel
model=[];
model.id=[];
model.name=[];
model.annotation=[];
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
model.rxnComps=[];
model.grRules={};
model.rxnGeneMat=[];
model.subSystems={};
model.eccodes={};
model.rxnMiriams={};
model.rxnNotes={};
model.rxnReferences={};
model.rxnConfidenceScores=[];
model.genes={};
model.geneComps=[];
model.geneMiriams={};
model.geneShortNames={};
model.proteins={};
model.metNames={};
model.metComps=[];
model.inchis={};
model.metFormulas={};
model.metMiriams={};
model.metCharges=[];
model.unconstrained=[];

%Load the model using libSBML
[modelSBML,errorMsg] = TranslateSBML_RAVEN(fileName,0,0,[1 1]);
if exist('tempFile','var')
    delete(tempFile)
end

if isempty(modelSBML)
    EM=['There is a problem with the SBML file. Try using the SBML Validator at http://sbml.org/Facilities/Validator.\nlibSBML reports: ', errorMsg.message];
    dispEM(EM);
end

%Retrieve compartment names and IDs
compartmentNames=cell(numel(modelSBML.compartment),1);
compartmentIDs=cell(numel(modelSBML.compartment),1);
compartmentOutside=cell(numel(modelSBML.compartment),1);
compartmentMiriams=cell(numel(modelSBML.compartment),1);

if isfield(modelSBML.compartment,'sboTerm') && numel(unique([modelSBML.compartment.sboTerm])) == 1
    %If all the SBO terms are identical, don't add them to compMiriams
    modelSBML.compartment = rmfield(modelSBML.compartment,'sboTerm');
end

for i=1:numel(modelSBML.compartment)
    compartmentNames{i}=modelSBML.compartment(i).name;
    compartmentIDs{i}=modelSBML.compartment(i).id;
    if isfield(modelSBML.compartment(i),'outside')
        if ~isempty(modelSBML.compartment(i).outside)
            compartmentOutside{i}=modelSBML.compartment(i).outside;
        else
            compartmentOutside{i}='';
        end
    else
        compartmentOutside{i}=[];
    end

    if isfield(modelSBML.compartment(i),'annotation')
        compartmentMiriams{i}=parseMiriam(modelSBML.compartment(i).annotation);
    else
        compartmentMiriams{i}=[];
    end

    if isfield(modelSBML.compartment(i),'sboTerm') && ~(modelSBML.compartment(i).sboTerm==-1)
        compartmentMiriams{i} = addSBOtoMiriam(compartmentMiriams{i},modelSBML.compartment(i).sboTerm);
    end
end

%If there are no compartment names then use compartment id as name
if all(cellfun(@isempty,compartmentNames))
    compartmentNames=compartmentIDs;
end

%Retrieve info on metabolites, genes, complexes
metaboliteNames={};
metaboliteIDs={};
metaboliteCompartments={};
metaboliteUnconstrained=[];
metaboliteFormula={};
metaboliteInChI={};
metaboliteMiriams={};
metaboliteCharges=[];

geneNames={};
geneIDs={};
geneMiriams={};
geneShortNames={};
proteins={};
geneCompartments={};
complexIDs={};
complexNames={};

%If the file is not a COBRA Toolbox model. According to the format
%specified in the yeast consensus model both metabolites and genes are a
%type of 'species'. The metabolites have names starting with 'M_' and genes
%with 'E_'
geneSBOs = [];
metSBOs = [];
%Regex of compartment names, later to be used to remove from metabolite
%names if present as suffix.
regexCompNames = ['\s?\[((' strjoin({modelSBML.compartment.name},')|(') '))\]$'];
for i=1:numel(modelSBML.species)
    if length(modelSBML.species(i).id)>=2 && strcmpi(modelSBML.species(i).id(1:2),'E_')
        geneNames{numel(geneNames)+1,1}=modelSBML.species(i).name;

        %The "E_" is included in the ID. This is because it's only used
        %internally in this file and it makes the matching a little
        %smoother
        geneIDs{numel(geneIDs)+1,1}=modelSBML.species(i).id;
        geneCompartments{numel(geneCompartments)+1,1}=modelSBML.species(i).compartment;

        %Get Miriam structure
        if isfield(modelSBML.species(i),'annotation')
            %Get Miriam info
            geneMiriam=parseMiriam(modelSBML.species(i).annotation);
            geneMiriams{numel(geneMiriams)+1,1}=geneMiriam;
        else
            geneMiriams{numel(geneMiriams)+1,1}=[];
        end

        %Protein short names (for example ERG10) are saved as SHORT
        %NAME: NAME in the notes-section of metabolites for SBML Level
        %2 and as PROTEIN_ASSOCIATION for each reaction in SBML Level 2
        %COBRA Toolbox format. For now only the SHORT NAME is loaded
        %and no mapping takes place
        if isfield(modelSBML.species(i),'notes')
            geneShortNames{numel(geneShortNames)+1,1}=parseNote(modelSBML.species(i).notes,'SHORT NAME');
        else
            geneShortNames{numel(geneShortNames)+1,1}='';
        end

        %Get SBO term
        if isfield(modelSBML.species(i),'sboTerm') && ~(modelSBML.species(i).sboTerm==-1)
            geneSBOs(end+1,1) = modelSBML.species(i).sboTerm;
        end
    elseif length(modelSBML.species(i).id)>=2 && strcmpi(modelSBML.species(i).id(1:3),'Cx_')
        %If it's a complex keep the ID and name
        complexIDs=[complexIDs;modelSBML.species(i).id];
        complexNames=[complexNames;modelSBML.species(i).name];
    else
        %If it is not gene or complex, then it must be a metabolite
        metaboliteNames{numel(metaboliteNames)+1,1}=modelSBML.species(i).name;
        metaboliteIDs{numel(metaboliteIDs)+1,1}=modelSBML.species(i).id;
        metaboliteCompartments{numel(metaboliteCompartments)+1,1}=modelSBML.species(i).compartment;
        metaboliteUnconstrained(numel(metaboliteUnconstrained)+1,1)=modelSBML.species(i).boundaryCondition;

        %For each metabolite retrieve the formula and the InChI code if
        %available First add the InChI code and the formula from the
        %InChI. This allows for overwriting the formula by setting the
        %actual formula field
        if ~isempty(modelSBML.species(i).annotation)
            %Get the formula if available
            startString='>InChI=';
            endString='</in:inchi>';
            formStart=strfind(modelSBML.species(i).annotation,startString);
            if isempty(formStart)
                startString='InChI=';
                endString='"/>';
            end
            formStart=strfind(modelSBML.species(i).annotation,startString);
            if ~isempty(formStart)
                formEnd=strfind(modelSBML.species(i).annotation,endString);
                formEndIndex=find(formEnd>formStart, 1 );
                formula=modelSBML.species(i).annotation(formStart+numel(startString):formEnd(formEndIndex)-1);
                metaboliteInChI{numel(metaboliteInChI)+1,1}=formula;

                %The composition is most often present between the
                %first and second "/" in the model. In some simple
                %molecules, such as salts, there is no second "/". The
                %formula is then assumed to be to the end of the string
                compositionIndexes=strfind(formula,'/');
                if numel(compositionIndexes)>1
                    metaboliteFormula{numel(metaboliteFormula)+1,1}=...
                        formula(compositionIndexes(1)+1:compositionIndexes(2)-1);
                else
                    if numel(compositionIndexes)==1
                        %Probably a simple molecule which can have only
                        %one conformation
                        metaboliteFormula{numel(metaboliteFormula)+1,1}=...
                            formula(compositionIndexes(1)+1:numel(formula));
                    else
                        metaboliteFormula{numel(metaboliteFormula)+1,1}='';
                    end
                end
            elseif isfield(modelSBML.species(i),'fbc_chemicalFormula')
                metaboliteInChI{numel(metaboliteInChI)+1,1}='';
                if ~isempty(modelSBML.species(i).fbc_chemicalFormula)
                    %Cannot extract InChi from formula, so remains
                    %empty
                    metaboliteFormula{numel(metaboliteFormula)+1,1}=modelSBML.species(i).fbc_chemicalFormula;
                else
                    metaboliteFormula{numel(metaboliteFormula)+1,1}='';
                end
            else
                metaboliteInChI{numel(metaboliteInChI)+1,1}='';
                metaboliteFormula{numel(metaboliteFormula)+1,1}='';
            end

            %Get Miriam info
            metMiriam=parseMiriam(modelSBML.species(i).annotation);
            metaboliteMiriams{numel(metaboliteMiriams)+1,1}=metMiriam;
        else
            metaboliteInChI{numel(metaboliteInChI)+1,1}='';
            if isfield(modelSBML.species(i),'notes')
                metaboliteFormula{numel(metaboliteFormula)+1,1}=parseNote(modelSBML.species(i).notes,'FORMULA');
            else
                metaboliteFormula{numel(metaboliteFormula)+1,1}='';
            end
            metaboliteMiriams{numel(metaboliteMiriams)+1,1}=[];
        end
        if ~isempty(modelSBML.species(i).notes)
            if ~isfield(modelSBML.species(i),'annotation')
                metaboliteFormula{numel(metaboliteFormula)+1,1}=parseNote(modelSBML.species(i).notes,'FORMULA');
            end
        elseif ~isfield(modelSBML.species(i),'annotation')
            metaboliteFormula{numel(metaboliteFormula)+1,1}='';
        end
        %Get SBO term
        if isfield(modelSBML.species(i),'sboTerm') && ~(modelSBML.species(i).sboTerm==-1)
            metSBOs(end+1,1) = modelSBML.species(i).sboTerm;
        end
    end

    %The following lines are executed regardless isSBML2COBRA setting
    if isempty(modelSBML.species(i).id) || ~strcmpi(modelSBML.species(i).id(1:2),'E_')
        if isempty(modelSBML.species(i).id) || ~strcmpi(modelSBML.species(i).id(1:3),'Cx_')
            %Remove trailing [compartment] from metabolite name if present
            metaboliteNames{end,1}=regexprep(metaboliteNames{end,1},regexCompNames,'');
            metaboliteNames{end,1}=metaboliteNames{end,1};
            if isfield(modelSBML.species(i),'fbc_charge')
                if ~isempty(modelSBML.species(i).fbc_charge) && modelSBML.species(i).isSetfbc_charge
                    metaboliteCharges(numel(metaboliteCharges)+1,1)=double(modelSBML.species(i).fbc_charge);
                else
                    if isfield(modelSBML.species(i),'notes')
                        if strfind(modelSBML.species(i).notes,'CHARGE')
                            metaboliteCharges(numel(metaboliteCharges)+1,1)=str2double(parseNote(modelSBML.species(i).notes,'CHARGE'));
                        else
                            metaboliteCharges(numel(metaboliteCharges)+1,1)=NaN;
                        end
                    else
                        metaboliteCharges(numel(metaboliteCharges)+1,1)=NaN;
                    end
                end
            elseif isfield(modelSBML.species(i),'notes')
                if strfind(modelSBML.species(i).notes,'CHARGE')
                    metaboliteCharges(numel(metaboliteCharges)+1,1)=str2double(parseNote(modelSBML.species(i).notes,'CHARGE'));
                else
                    metaboliteCharges(numel(metaboliteCharges)+1,1)=NaN;
                end
            else
                metaboliteCharges(numel(metaboliteCharges)+1,1)=NaN;
            end
            %Additional information from FBC format Chemical formula
            if isfield(modelSBML.species(i),'fbc_chemicalFormula')
                if ~isempty(modelSBML.species(i).fbc_chemicalFormula)
                    metaboliteFormula{numel(metaboliteFormula),1}=modelSBML.species(i).fbc_chemicalFormula;
                end
            end
        end
    end
end

%Add SBO terms to gene and metabolite miriam fields
if numel(unique(geneSBOs)) > 1  % don't add if they're all identical
    for i = 1:numel(geneNames)
        geneMiriams{i} = addSBOtoMiriam(geneMiriams{i},geneSBOs(i));
    end
end
if numel(unique(metSBOs)) > 1
    for i = 1:numel(metaboliteNames)
        metaboliteMiriams{i} = addSBOtoMiriam(metaboliteMiriams{i},metSBOs(i));
    end
end

%Retrieve info on reactions
reactionNames=cell(numel(modelSBML.reaction),1);
reactionIDs=cell(numel(modelSBML.reaction),1);
subsystems=cell(numel(modelSBML.reaction),1);
eccodes=cell(numel(modelSBML.reaction),1);
eccodes(:,:)=cellstr('');
rxnconfidencescores=NaN(numel(modelSBML.reaction),1);
rxnreferences=cell(numel(modelSBML.reaction),1);
rxnreferences(:,:)=cellstr('');
rxnnotes=cell(numel(modelSBML.reaction),1);
rxnnotes(:,:)=cellstr('');
grRules=cell(numel(modelSBML.reaction),1);
grRules(:,:)=cellstr('');
grRulesFromModifier=grRules;
rxnComps=zeros(numel(modelSBML.reaction),1);
rxnMiriams=cell(numel(modelSBML.reaction),1);
reactionReversibility=zeros(numel(modelSBML.reaction),1);
reactionUB=zeros(numel(modelSBML.reaction),1);
reactionLB=zeros(numel(modelSBML.reaction),1);
reactionObjective=zeros(numel(modelSBML.reaction),1);

%Construct the stoichiometric matrix while the reaction info is read
S=zeros(numel(metaboliteIDs),numel(modelSBML.reaction));

counter=0;
%If FBC, then bounds have parameter ids defined for the whole model
if isfield(modelSBML,'parameter')
    parameter.name=cell(numel(modelSBML.parameter),1);
    parameter.name={modelSBML.parameter(:).id}';
    parameter.value={modelSBML.parameter(:).value}';
end

if isfield(modelSBML.reaction,'sboTerm') && numel(unique([modelSBML.reaction.sboTerm])) == 1
    %If all the SBO terms are identical, don't add them to rxnMiriams
    modelSBML.reaction = rmfield(modelSBML.reaction,'sboTerm');
end

for i=1:numel(modelSBML.reaction)

    %Check that the reaction doesn't produce a complex and nothing else. If
    %so, then jump to the next reaction. This is because I get the genes
    %for complexes from the names and not from the reactions that create
    %them. This only applies to the non-COBRA format
    if numel(modelSBML.reaction(i).product)==1
        if length(modelSBML.reaction(i).product(1).species)>=3
            if strcmp(modelSBML.reaction(i).product(1).species(1:3),'Cx_')==true
                continue;
            end
        end
    end

    %It didn't look like a gene complex-forming reaction
    counter=counter+1;

    reactionNames{counter}=modelSBML.reaction(i).name;

    reactionIDs{counter}=modelSBML.reaction(i).id;
    reactionReversibility(counter)=modelSBML.reaction(i).reversible;

    %If model is FBC, first get parameter of bound and then replace it with
    %the correct value. Probably faster with replace(), but this was only
    %introduced in Matlab R2016b
    if isfield(modelSBML.reaction(i),'fbc_lowerFluxBound')
        lb=modelSBML.reaction(i).fbc_lowerFluxBound;
        ub=modelSBML.reaction(i).fbc_upperFluxBound;
        for n=1:numel(parameter.value)
            lb=regexprep(lb,parameter.name(n),num2str(parameter.value{n}));
            ub=regexprep(ub,parameter.name(n),num2str(parameter.value{n}));
        end
        if isempty(lb)
            lb='-Inf';
        end
        if isempty(ub)
            ub='Inf';
        end
        reactionLB(counter)=str2num(lb);
        reactionUB(counter)=str2num(ub);
        %The order of these parameters should not be hard coded
    elseif isfield(modelSBML.reaction(i).kineticLaw,'parameter')
        reactionLB(counter)=modelSBML.reaction(i).kineticLaw.parameter(1).value;
        reactionUB(counter)=modelSBML.reaction(i).kineticLaw.parameter(2).value;
        reactionObjective(counter)=modelSBML.reaction(i).kineticLaw.parameter(3).value;
    else
        if reactionReversibility(counter)==true
            reactionLB(counter)=-inf;
        else
            reactionLB(counter)=0;
        end
        reactionUB(counter)=inf;
        reactionObjective(counter)=0;
    end

    %Find the associated gene if available
    %If FBC, get gene association data from corresponding fields
    if isfield(modelSBML.reaction(i),'fbc_geneProductAssociation')
        if ~isempty(modelSBML.reaction(i).fbc_geneProductAssociation) && ~isempty(modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association)
            grRules{counter}=modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association.fbc_association;
        end
    elseif isfield(modelSBML.reaction(i),'notes')
        %This section was previously executed only if isSBML2COBRA is true. Now
        %it will be executed, if 'GENE_ASSOCIATION' is found in
        %modelSBML.reaction(i).notes
        if strfind(modelSBML.reaction(i).notes,'GENE_ASSOCIATION')
            geneAssociation=parseNote(modelSBML.reaction(i).notes,'GENE_ASSOCIATION');
        elseif strfind(modelSBML.reaction(i).notes,'GENE ASSOCIATION')
            geneAssociation=parseNote(modelSBML.reaction(i).notes,'GENE ASSOCIATION');
        else
            geneAssociation='';
        end
        if ~isempty(geneAssociation)
            %This adds the grRules. The gene list and rxnGeneMat are created
            %later
            grRules{counter}=geneAssociation;
        end
    end
    if isempty(grRules{counter}) && ~isempty(modelSBML.reaction(i).modifier)
        rules='';
        for j=1:numel(modelSBML.reaction(i).modifier)
            modifier=modelSBML.reaction(i).modifier(j).species;
            if ~isempty(modifier)
                if strcmpi(modifier(1:2),'E_')
                    index=find(strcmp(modifier,geneIDs));
                    %This should be unique and in the geneIDs list,
                    %otherwise something is wrong
                    if numel(index)~=1
                        EM=['Could not get the gene association data from reaction ' reactionIDs{i}];
                        dispEM(EM);
                    end
                    if ~isempty(rules)
                        rules=[rules ' or (' geneNames{index} ')'];
                    else
                        rules=['(' geneNames{index} ')'];
                    end
                elseif strcmp(modifier(1:2),'s_')
                    index=find(strcmp(modifier,metaboliteIDs));
                    %This should be unique and in the geneIDs list,
                    %otherwise something is wrong
                    if numel(index)~=1
                        EM=['Could not get the gene association data from reaction ' reactionIDs{i}];
                        dispEM(EM);
                    end
                    if ~isempty(rules)
                        rules=[rules ' or (' metaboliteIDs{index} ')'];
                    else
                        rules=['(' metaboliteIDs{index} ')'];
                    end
                else
                    %It seems to be a complex. Add the corresponding
                    %genes from the name of the complex (not the
                    %reaction that creates it)
                    index=find(strcmp(modifier,complexIDs));
                    if numel(index)==1
                        if ~isempty(rules)
                            rules=[rules ' or (' strrep(complexNames{index},':',' and ') ')'];
                        else
                            rules=['(' strrep(complexNames{index},':',' and ') ')'];
                        end
                    else
                        %Could not find a complex
                        EM=['Could not get the gene association data from reaction ' reactionIDs{i}];
                        dispEM(EM);
                    end
                end
            end
        end
        grRules{counter}=rules;
        grRulesFromModifier{counter}=rules;%Backup copy for grRules, useful to parse Yeast 7.6
    end

    %Add reaction compartment
    if isfield(modelSBML.reaction(i),'compartment')
        if ~isempty(modelSBML.reaction(i).compartment)
            rxnComp=modelSBML.reaction(i).compartment;
        else
            rxnComp='';
        end
    elseif isfield(modelSBML.reaction(i),'notes')
        rxnComp=parseNote(modelSBML.reaction(i).notes,'COMPARTMENT');
    end
    if ~isempty(rxnComp)
        %Find it in the compartment list
        [~, J]=ismember(rxnComp,compartmentIDs);
        rxnComps(counter)=J;
    end


    miriamStruct=parseMiriam(modelSBML.reaction(i).annotation);
    rxnMiriams{counter}=miriamStruct;
    if isfield(modelSBML.reaction(i),'notes')
        subsystems{counter,1}=cellstr(parseNote(modelSBML.reaction(i).notes,'SUBSYSTEM'));
        subsystems{counter,1}(cellfun('isempty',subsystems{counter,1})) = [];
        if strfind(modelSBML.reaction(i).notes,'Confidence Level')
            confScore = parseNote(modelSBML.reaction(i).notes,'Confidence Level');
            if isempty(confScore)
                confScore = 0;
            end
            rxnconfidencescores(counter)=str2double(confScore);
        end
        rxnreferences{counter,1}=parseNote(modelSBML.reaction(i).notes,'AUTHORS');
        rxnnotes{counter,1}=parseNote(modelSBML.reaction(i).notes,'NOTES');
    end

    %Get SBO terms
    if isfield(modelSBML.reaction(i),'sboTerm') && ~(modelSBML.reaction(i).sboTerm==-1)
        rxnMiriams{counter} = addSBOtoMiriam(rxnMiriams{counter}, modelSBML.reaction(i).sboTerm);
    end

    %Get ec-codes
    eccode='';
    if ~isempty(modelSBML.reaction(i).annotation)
        if strfind(modelSBML.reaction(i).annotation,'urn:miriam:ec-code')
            eccode=parseAnnotation(modelSBML.reaction(i).annotation,'urn:miriam:',':','ec-code');
        elseif strfind(modelSBML.reaction(i).annotation,'http://identifiers.org/ec-code')
            eccode=parseAnnotation(modelSBML.reaction(i).annotation,'http://identifiers.org/','/','ec-code');
        elseif strfind(modelSBML.reaction(i).annotation,'https://identifiers.org/ec-code')
            eccode=parseAnnotation(modelSBML.reaction(i).annotation,'https://identifiers.org/','/','ec-code');
        end
    elseif isfield(modelSBML.reaction(i),'notes')
        if strfind(modelSBML.reaction(i).notes,'EC Number')
            eccode=[eccode parseNote(modelSBML.reaction(i).notes,'EC Number')];
        elseif strfind(modelSBML.reaction(i).notes,'PROTEIN_CLASS')
            eccode=[eccode parseNote(modelSBML.reaction(i).notes,'PROTEIN_CLASS')];
        end
    end
    eccodes{counter}=eccode;

    %Add all reactants
    for j=1:numel(modelSBML.reaction(i).reactant)
        %Get the index of the metabolite in metaboliteIDs. External
        %metabolites will be removed at a later stage
        metIndex=find(strcmp(modelSBML.reaction(i).reactant(j).species,metaboliteIDs),1);
        if isempty(metIndex)
            EM=['Could not find metabolite ' modelSBML.reaction(i).reactant(j).species ' in reaction ' reactionIDs{counter}];
            dispEM(EM);
        end
        S(metIndex,counter)=S(metIndex,counter)+modelSBML.reaction(i).reactant(j).stoichiometry*-1;
    end

    %Add all products
    for j=1:numel(modelSBML.reaction(i).product)
        %Get the index of the metabolite in metaboliteIDs.
        metIndex=find(strcmp(modelSBML.reaction(i).product(j).species,metaboliteIDs),1);
        if isempty(metIndex)
            EM=['Could not find metabolite ' modelSBML.reaction(i).product(j).species ' in reaction ' reactionIDs{counter}];
            dispEM(EM);
        end
        S(metIndex,counter)=S(metIndex,counter)+modelSBML.reaction(i).product(j).stoichiometry;
    end
end

%if FBC, objective function is separately defined. Multiple objective
%functions can be defined, one is set as active
if isfield(modelSBML, 'fbc_activeObjective')
    obj=modelSBML.fbc_activeObjective;
    for i=1:numel(modelSBML.fbc_objective)
        if strcmp(obj,modelSBML.fbc_objective(i).fbc_id)
            if ~isempty(modelSBML.fbc_objective(i).fbc_fluxObjective)
                rxn=modelSBML.fbc_objective(i).fbc_fluxObjective.fbc_reaction;
                idx=find(ismember(reactionIDs,rxn));
                reactionObjective(idx)=modelSBML.fbc_objective(i).fbc_fluxObjective.fbc_coefficient;
            end
        end
    end
end

%subSystems can be stored as groups instead of in annotations
if isfield(modelSBML,'groups_group')
    for i=1:numel(modelSBML.groups_group)
        groupreactions={modelSBML.groups_group(i).groups_member(:).groups_idRef};
        [~, idx] = ismember(groupreactions, reactionIDs);
        if any(idx)
            for j=1:numel(idx)
                if isempty(subsystems{idx(j)}) % First subsystem
                    subsystems{idx(j)} = {modelSBML.groups_group(i).groups_name};
                else % Consecutive subsystems: concatenate
                    subsystems{idx(j)} = horzcat(subsystems{idx(j)}, modelSBML.groups_group(i).groups_name);
                end
            end
        end
    end
end

%Shrink the structures if complex-forming reactions had to be skipped
reactionNames=reactionNames(1:counter);
reactionIDs=reactionIDs(1:counter);
subsystems=subsystems(1:counter);
eccodes=eccodes(1:counter);
rxnconfidencescores=rxnconfidencescores(1:counter);
rxnreferences=rxnreferences(1:counter);
rxnnotes=rxnnotes(1:counter);
grRules=grRules(1:counter);
rxnMiriams=rxnMiriams(1:counter);
reactionReversibility=reactionReversibility(1:counter);
reactionUB=reactionUB(1:counter);
reactionLB=reactionLB(1:counter);
reactionObjective=reactionObjective(1:counter);
S=S(:,1:counter);

model.name=modelSBML.name;
model.id=modelSBML.id;
model.rxns=reactionIDs;
model.mets=metaboliteIDs;
model.S=sparse(S);
model.lb=reactionLB;
model.ub=reactionUB;
model.rev=reactionReversibility;
model.c=reactionObjective;
model.b=zeros(numel(metaboliteIDs),1);
model.comps=compartmentIDs;
model.compNames=compartmentNames;
model.rxnConfidenceScores=rxnconfidencescores;
model.rxnReferences=rxnreferences;
model.rxnNotes=rxnnotes;

%Load annotation if available. If there are several authors, only the first
%author credentials are imported
if isfield(modelSBML,'annotation')
    endString='</';
    I=strfind(modelSBML.annotation,endString);
    J=strfind(modelSBML.annotation,'<vCard:Family>');
    if any(J)
        model.annotation.familyName=modelSBML.annotation(J(1)+14:I(find(I>J(1),1))-1);
    end
    J=strfind(modelSBML.annotation,'<vCard:Given>');
    if any(J)
        model.annotation.givenName=modelSBML.annotation(J(1)+13:I(find(I>J(1),1))-1);
    end
    J=strfind(modelSBML.annotation,'<vCard:EMAIL>');
    if any(J)
        model.annotation.email=modelSBML.annotation(J(1)+13:I(find(I>J(1),1))-1);
    end
    J=strfind(modelSBML.annotation,'<vCard:Orgname>');
    if any(J)
        model.annotation.organization=modelSBML.annotation(J(1)+15:I(find(I>J(1),1))-1);
    end
    endString='"/>';
    I=strfind(modelSBML.annotation,endString);
    if strfind(modelSBML.annotation,'"urn:miriam:')
        J=strfind(modelSBML.annotation,'"urn:miriam:');
        if any(J)
            model.annotation.taxonomy=modelSBML.annotation(J+12:I(find(I>J,1))-1);
        end
    else
        J=strfind(modelSBML.annotation,'"http://identifiers.org/');
        if any(J)
            model.annotation.taxonomy=modelSBML.annotation(J+24:I(find(I>J,1))-1);
        else
            J=strfind(modelSBML.annotation,'"https://identifiers.org/');
            if any(J)
                model.annotation.taxonomy=modelSBML.annotation(J+25:I(find(I>J,1))-1);
            end
        end
    end
end
if isfield(modelSBML,'notes')
    startString=strfind(modelSBML.notes,'xhtml">');
    endString=strfind(modelSBML.notes,'</body>');
    if any(startString) && any(endString)
        model.annotation.note=modelSBML.notes(startString+7:endString-1);
        model.annotation.note=regexprep(model.annotation.note,'<p>|</p>','');
        model.annotation.note=strtrim(model.annotation.note);
        if regexp(model.annotation.note,'This file was generated using the exportModel function in RAVEN Toolbox \d\.\d and OutputSBML in libSBML')
            model.annotation=rmfield(model.annotation,'note'); % Default note added when running exportModel
        end
    end
end

if any(~cellfun(@isempty,compartmentOutside))
    model.compOutside=compartmentOutside;
end

model.rxnNames=reactionNames;
model.metNames=metaboliteNames;

%Match the compartments for metabolites
[~, J]=ismember(metaboliteCompartments,model.comps);
model.metComps=J;

%If any genes have been loaded (only for the new format)
if ~isempty(geneNames)
    %In some rare cases geneNames may not necessarily be used in grRules.
    %That is true for Yeast 7.6. It's therefore important to change gene
    %systematic names to geneIDs in sophisticated way. Gene systematic
    %names are not unique, since exactly the same name may be in different
    %compartments
    if all(cellfun(@isempty,strfind(grRules,geneNames{1})))
        geneShortNames=geneNames;
        %geneShortNames contain compartments as well, so these are removed
        geneShortNames=regexprep(geneShortNames,' \[.+$','');
        %grRules obtained from modifier fields contain geneNames. These are
        %changed into geneIDs. grRulesFromModifier is a good way to have
        %geneIDs and rxns association when it's important to resolve
        %systematic name ambiguities
        grRulesFromModifier=regexprep(regexprep(grRulesFromModifier,'\[|\]','_'),regexprep(geneNames,'\[|\]','_'),geneIDs);
        grRules=regexprep(regexprep(grRules,'\[|\]','_'),regexprep(geneNames,'\[|\]','_'),geneIDs);

        %Yeast 7.6 contains several metabolites, which were used in gene
        %associations. For that reason, the list of species ID is created
        %and we then check whether any of them have kegg.genes annotation
        %thereby obtaining systematic gene names
        geneShortNames=vertcat(geneShortNames,metaboliteNames);
        geneIDs=vertcat(geneIDs,metaboliteIDs);
        geneSystNames=extractMiriam(vertcat(geneMiriams,metaboliteMiriams),'kegg.genes');
        geneCompartments=vertcat(geneCompartments,metaboliteCompartments);
        geneMiriams=vertcat(geneMiriams,metaboliteMiriams);

        %Now we retain information for only these entries, which have
        %kegg.genes annotation
        geneShortNames=geneShortNames(~cellfun('isempty',geneSystNames));
        geneIDs=geneIDs(~cellfun('isempty',geneSystNames));
        geneSystNames=geneSystNames(~cellfun('isempty',geneSystNames));
        geneCompartments=geneCompartments(~cellfun('isempty',geneSystNames));
        geneMiriams=geneMiriams(~cellfun('isempty',geneSystNames));
        %Now we reorder geneIDs and geneSystNames by geneSystNames string
        %length
        geneNames=geneIDs;%Backuping geneIDs, since we need unsorted order for later
        [~, Indx] = sort(cellfun('size', geneSystNames, 2), 'descend');
        geneIDs = geneIDs(Indx);
        geneSystNames = geneSystNames(Indx);
        for i=1:numel(geneSystNames)
            for j=1:numel(grRules)
                if strfind(grRules{j},geneSystNames{i})
                    if ~isempty(grRules{j})
                        if sum(ismember(geneSystNames,geneSystNames{i}))==1
                            grRules{j}=regexprep(grRules{j},geneSystNames{i},geneIDs{i});
                        elseif sum(ismember(geneSystNames,geneSystNames{i}))>1
                            counter=0;
                            ovrlpIDs=geneIDs(ismember(geneSystNames,geneSystNames{i}));
                            for k=1:numel(ovrlpIDs)
                                if strfind(grRulesFromModifier{j},ovrlpIDs{k})
                                    counter=counter+1;
                                    grRules{j}=regexprep(grRules{j},geneSystNames{i},ovrlpIDs{k});
                                end
                                if counter>1
                                    EM=['Gene association is ambiguous for reaction ' modelSBML.reaction(j).id];
                                    dispEM(EM);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    model.genes=geneNames;
    model.grRules=grRules;
    [grRules,rxnGeneMat] = standardizeGrRules(model,true);
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;

    %Match the compartments for genes
    [~, J]=ismember(geneCompartments,model.comps);
    model.geneComps=J;
else
    if ~all(cellfun(@isempty,grRules))
        %If fbc_geneProduct exists, follow the specified gene order, such
        %that matching geneShortNames in function below will work
        if isfield(modelSBML,'fbc_geneProduct')
            genes={modelSBML.fbc_geneProduct.fbc_id};

            %Get gene Miriams if they were not retrieved above (this occurs
            %when genes are stored as fbc_geneProduct instead of species)
            if isempty(geneMiriams)
                geneMiriams = cell(numel(genes),1);
                if isfield(modelSBML.fbc_geneProduct,'sboTerm') && numel(unique([modelSBML.fbc_geneProduct.sboTerm])) == 1
                    %If all the SBO terms are identical, don't add them to geneMiriams
                    modelSBML.fbc_geneProduct = rmfield(modelSBML.fbc_geneProduct,'sboTerm');
                end
                for i = 1:numel(genes)
                    geneMiriams{i}=parseMiriam(modelSBML.fbc_geneProduct(i).annotation);
                    if isfield(modelSBML.fbc_geneProduct(i),'sboTerm') && ~(modelSBML.fbc_geneProduct(i).sboTerm==-1)
                        geneMiriams{i} = addSBOtoMiriam(geneMiriams{i},modelSBML.fbc_geneProduct(i).sboTerm);
                    end
                end
            end
            proteins={modelSBML.fbc_geneProduct.fbc_name};
        else
            genes=getGeneList(grRules);
        end
        model.genes=genes;
        model.grRules=grRules;
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        model.grRules = grRules;
        model.rxnGeneMat = rxnGeneMat;
    end
end

if all(cellfun(@isempty,geneShortNames))
    if isfield(modelSBML,'fbc_geneProduct')
        for i=1:numel(genes)
            if ~isempty(modelSBML.fbc_geneProduct(i).fbc_label)
                geneShortNames{i,1}=modelSBML.fbc_geneProduct(i).fbc_label;
            elseif ~isempty(modelSBML.fbc_geneProduct(i).fbc_name)
                geneShortNames{i,1}=modelSBML.fbc_geneProduct(i).fbc_name;
            else
                geneShortNames{i,1}='';
            end
        end
    end
end

%If any InChIs have been loaded
if any(~cellfun(@isempty,metaboliteInChI))
    model.inchis=metaboliteInChI;
end

%If any formulas have been loaded
if any(~cellfun(@isempty,metaboliteFormula))
    model.metFormulas=metaboliteFormula;
end

%If any charges have been loaded
if ~isempty(metaboliteCharges)
    model.metCharges=metaboliteCharges;
end

%If any gene short names have been loaded
if any(~cellfun(@isempty,geneShortNames))
    model.geneShortNames=geneShortNames;
end

%If any Miriam strings for compartments have been loaded
if any(~cellfun(@isempty,compartmentMiriams))
    model.compMiriams=compartmentMiriams;
end

%If any Miriam strings for metabolites have been loaded
if any(~cellfun(@isempty,metaboliteMiriams))
    model.metMiriams=metaboliteMiriams;
end

%If any subsystems have been loaded
if any(~cellfun(@isempty,subsystems))
    model.subSystems=subsystems;
end
if any(rxnComps)
    if all(rxnComps)
        model.rxnComps=rxnComps;
    else
        if supressWarnings==false
            EM='The compartments for the following reactions could not be matched. Ignoring reaction compartment information';
            dispEM(EM,false,model.rxns(rxnComps==0));
        end
    end
end

%If any ec-codes have been loaded
if any(~cellfun(@isempty,eccodes))
    model.eccodes=eccodes;
end

%If any Miriam strings for reactions have been loaded
if any(~cellfun(@isempty,rxnMiriams))
    model.rxnMiriams=rxnMiriams;
end

%If any Miriam strings for genes have been loaded
if any(~cellfun(@isempty,geneMiriams))
    model.geneMiriams=geneMiriams;
end

%If any protein strings have been loaded
if any(~cellfun(@isempty,proteins))
    proteins = reshape(proteins,[],1);
    model.proteins=proteins;
end

model.unconstrained=metaboliteUnconstrained;

%Convert SBML IDs back into their original strings. Here we are using part
%from convertSBMLID, originating from the COBRA Toolbox
model.rxns=regexprep(model.rxns,'__([0-9]+)__','${char(str2num($1))}');
model.mets=regexprep(model.mets,'__([0-9]+)__','${char(str2num($1))}');
model.comps=regexprep(model.comps,'__([0-9]+)__','${char(str2num($1))}');
model.grRules=regexprep(model.grRules,'__([0-9]+)__','${char(str2num($1))}');
model.genes=regexprep(model.genes,'__([0-9]+)__','${char(str2num($1))}');
model.id=regexprep(model.id,'__([0-9]+)__','${char(str2num($1))}');


hasCOBRAids(1) = all(startsWith(model.rxns,'R_'));
hasCOBRAids(2) = all(startsWith(model.mets,'M_'));
hasCOBRAids(3) = all(startsWith(model.comps,'C_'));
hasCOBRAids(4) = all(startsWith(model.genes,'G_'));
hasCOBRAids(5) = all(startsWith(model.metNames,'M_'));
hasCOBRAids(6) = all(startsWith(model.rxnNames,'R_'));
hasCOBRAids(7) = all(startsWith(model.id,'M_'));

if COBRAstyle
    if hasCOBRAids(1)
        model.rxns  = cellfun(@(x) x(3:end), model.rxns,'un',0);
    end
    if hasCOBRAids(2)
        model.mets  = cellfun(@(x) x(3:end), model.mets,'un',0);
    end
    if hasCOBRAids(3)
        model.comps = cellfun(@(x) x(3:end), model.comps,'un',0);
    end
    if hasCOBRAids(4)
        model.genes = cellfun(@(x) x(3:end), model.genes,'un',0);
        model.grRules=regexprep(model.grRules,'^G_','');
        model.grRules=regexprep(model.grRules,'\(G_','(');
        model.grRules=regexprep(model.grRules,' G_',' ');
    end
    % Not strictly identifier, but also remove from metNames and rxnNames if present
    if hasCOBRAids(5)
        model.metNames  = cellfun(@(x) x(3:end), model.metNames,'un',0);
    end
    if hasCOBRAids(6)
        model.rxnNames  = cellfun(@(x) x(3:end), model.rxnNames,'un',0);
    end
    if hasCOBRAids(6)
        model.id        = model.id(3:end);
    end
elseif any(hasCOBRAids)
    hasCOBRAidsMsg = {'model.rxns (R_ prefix)',...
                      'model.mets (M_ prefix)',...
                      'model.comps (C_ prefix)',...
                      'model.genes (G_ prefix)',...
                      'model.metNames (M_ prefix)',...
                      'model.rxnNames (R_ prefix)',...
                      'model.id (M_ prefix)'};
    hasCOBRAidsMsg(~hasCOBRAids) = [];
    printOrange(['COBRA style identifier prefixes are found in: "' ...
        strjoin(hasCOBRAidsMsg,'"; "') '". Since RAVEN 2.10.0, identifiers ' ...
        'are by default imported "as-is". If you do prefer to remove these ' ...
        'identifier prefixes, run importModel with COBRAstyle as true. ' ...
        'Example:   importModel(''filename.xml'', [], true);\n']);
end

%Remove unused fields
if isempty(model.annotation)
    model=rmfield(model,'annotation');
end
if isempty(model.compOutside)
    model=rmfield(model,'compOutside');
end
if isempty(model.compMiriams)
    model=rmfield(model,'compMiriams');
end
if isempty(model.rxnComps)
    model=rmfield(model,'rxnComps');
end
if isempty(model.grRules)
    model=rmfield(model,'grRules');
end
if isempty(model.rxnGeneMat)
    model=rmfield(model,'rxnGeneMat');
end
if isempty(model.subSystems)
    model=rmfield(model,'subSystems');
else
    model.subSystems(cellfun(@isempty,subsystems))={{''}};
end
if isempty(model.eccodes)
    model=rmfield(model,'eccodes');
end
if isempty(model.rxnMiriams)
    model=rmfield(model,'rxnMiriams');
end
if cellfun(@isempty,model.rxnNotes)
    model=rmfield(model,'rxnNotes');
end
if cellfun(@isempty,model.rxnReferences)
    model=rmfield(model,'rxnReferences');
end
if isempty(model.rxnConfidenceScores) || all(isnan(model.rxnConfidenceScores))
    model=rmfield(model,'rxnConfidenceScores');
end
if isempty(model.genes)
    model=rmfield(model,'genes');
elseif isrow(model.genes)
    model.genes=transpose(model.genes);
end
if isempty(model.geneComps)
    model=rmfield(model,'geneComps');
end
if isempty(model.geneMiriams)
    model=rmfield(model,'geneMiriams');
end
if isempty(model.geneShortNames)
    model=rmfield(model,'geneShortNames');
end
if isempty(model.proteins)
    model=rmfield(model,'proteins');
end
if isempty(model.inchis)
    model=rmfield(model,'inchis');
end
if isempty(model.metFormulas)
    model=rmfield(model,'metFormulas');
end
if isempty(model.metMiriams)
    model=rmfield(model,'metMiriams');
end
if ~any(model.metCharges)
    model=rmfield(model,'metCharges');
end

%This just removes the grRules if no genes have been loaded
if ~isfield(model,'genes') && isfield(model,'grRules')
    model=rmfield(model,'grRules');
end

%Print warnings about bad structure
if supressWarnings==false
    checkModelStruct(model,false);
end

if removeExcMets==true
    model=simplifyModel(model);
end
end

function matchGenes=getGeneList(grRules)
%Constructs the list of unique genes from grRules

%Assumes that everything that isn't a paranthesis, " AND " or " or " is a
%gene name
genes=strrep(grRules,'(','');
genes=strrep(genes,')','');
genes=strrep(genes,' or ',' ');
genes=strrep(genes,' and ',' ');
genes=strrep(genes,' OR ',' ');
genes=strrep(genes,' AND ',' ');
genes=regexp(genes,' ','split');

allNames={};
for i=1:numel(genes)
    allNames=[allNames genes{i}];
end
matchGenes=unique(allNames)';

%Remove the empty element if present
if isempty(matchGenes{1})
    matchGenes(1)=[];
end
end

function fieldContent=parseNote(searchString,fieldName)
%The function obtains the particular information from 'notes' field, using
%fieldName as the dummy string

fieldContent='';

if strfind(searchString,fieldName)
    [~,targetString] = regexp(searchString,['<p>' fieldName '.*?</p>'],'tokens','match');
    targetString=regexprep(targetString,'<p>|</p>','');
    targetString=regexprep(targetString,[fieldName, ':'],'');
    for i=1:numel(targetString)
        fieldContent=[fieldContent ';' strtrim(targetString{1,i})];
    end
    fieldContent=regexprep(fieldContent,'^;|;$','');
else
    fieldContent='';
end
end

function fieldContent=parseAnnotation(searchString,startString,midString,fieldName)

fieldContent='';

%Removing whitespace characters from the ending strings, which may occur in
%several cases
searchString=regexprep(searchString,'" />','"/>');
[~,targetString] = regexp(searchString,['<rdf:li rdf:resource="' startString fieldName midString '.*?"/>'],'tokens','match');
targetString=regexprep(targetString,'<rdf:li rdf:resource="|"/>','');
targetString=regexprep(targetString,startString,'');
targetString=regexprep(targetString,[fieldName midString],'');

for i=1:numel(targetString)
    fieldContent=[fieldContent ';' strtrim(targetString{1,i})];
end

fieldContent=regexprep(fieldContent,'^;|;$','');
end

function miriamStruct=parseMiriam(searchString)
%Generates miriam structure from annotation field

%Finding whether miriams are written in the old or the new way
if strfind(searchString,'urn:miriam:')
    startString='urn:miriam:';
    midString=':';
elseif strfind(searchString,'http://identifiers.org/')
    startString='http://identifiers.org/';
    midString='/';
elseif strfind(searchString,'https://identifiers.org/')
    startString='https://identifiers.org/';
    midString='/';
else
    miriamStruct=[];
    return;
end

miriamStruct=[];

searchString=regexprep(searchString,'" />','"/>');
[~,targetString] = regexp(searchString,'<rdf:li rdf:resource=".*?"/>','tokens','match');
targetString=regexprep(targetString,'<rdf:li rdf:resource="|"/>','');
targetString=regexprep(targetString,startString,'');
targetString=regexprep(targetString,midString,'/','once');

counter=0;
for i=1:numel(targetString)
    if isempty(regexp(targetString{1,i},'inchi|ec-code', 'once'))
        counter=counter+1;
        miriamStruct.name{counter,1} = regexprep(targetString{1,i},'/.+','','once');
        miriamStruct.value{counter,1} = regexprep(targetString{1,i},[miriamStruct.name{counter,1} '/'],'','once');
        miriamStruct.name{counter,1} = regexprep(miriamStruct.name{counter,1},'^obo\.','');
    end
end
end

function miriam = addSBOtoMiriam(miriam,sboTerm)
%Appends SBO term to miriam structure

sboTerm = {['SBO:' sprintf('%07u',sboTerm)]};  % convert to proper format
if isempty(miriam)
    miriam.name = {'sbo'};
    miriam.value = sboTerm;
elseif any(strcmp('sbo',miriam.name))
    currSbo = strcmp('sbo',miriam.name);
    miriam.value(currSbo) = sboTerm;
else
    miriam.name(end+1) = {'sbo'};
    miriam.value(end+1) = sboTerm;
end
end
