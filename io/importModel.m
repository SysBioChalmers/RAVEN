function model=importModel(varargin)
% importModel  Import a constraint-based model from an SBML file.
%
% Name-Value Arguments
% --------------------
% fileName : char
%     a SBML file to import. A dialog window will open if no file name is
%     specified.
% removeExcMets : logical
%     if true, exchange metabolites are removed from the imported model
%     (default true).
% removePrefix : logical
%     true if identifier prefixes should be removed when loading the model:
%     G_ for genes, R_ for reactions, M_ for metabolites, and C_ for
%     compartments. These are only removed if all identifiers of a certain
%     type contain the prefix (default true).
% supressWarnings : logical
%     true if warnings regarding the model structure should be supressed
%     (default false).
%
% Returns
% -------
% model : struct
%     imported model structure with fields:
%
%     - id : model ID
%     - name : name of model contents
%     - annotation : additional information about model
%     - rxns : reaction ids
%     - mets : metabolite ids
%     - S : stoichiometric matrix
%     - lb : lower bounds
%     - ub : upper bounds
%     - rev : reversibility vector
%     - c : objective coefficients
%     - b : equality constraints for the metabolite equations
%     - comps : compartment ids
%     - compNames : compartment names
%     - compOutside : the id (as in comps) for the compartment surrounding
%       each of the compartments
%     - compMiriams : structure with MIRIAM information about the
%       compartments
%     - rxnNames : reaction description
%     - rxnComps : compartments for reactions
%     - grRules : reaction to gene rules in text form
%     - rxnGeneMat : reaction-to-gene mapping in sparse matrix form
%     - subSystems : subsystem name for each reaction
%     - eccodes : EC-codes for the reactions
%     - rxnMiriams : structure with MIRIAM information about the reactions
%     - rxnNotes : reaction notes
%     - rxnReferences : reaction references
%     - rxnConfidenceScores : reaction confidence scores
%     - genes : list of all genes
%     - geneComps : compartments for genes
%     - geneMiriams : structure with MIRIAM information about the genes
%     - geneShortNames : gene alternative names (e.g. ERG10)
%     - proteins : protein associated to each gene
%     - metNames : metabolite description
%     - metComps : compartments for metabolites
%     - inchis : InChI-codes for metabolites
%     - metFormulas : metabolite chemical formula
%     - metMiriams : structure with MIRIAM information about the metabolites
%     - metCharges : metabolite charge
%     - unconstrained : true if the metabolite is an exchange metabolite
%
% Examples
% --------
%     model = importModel(fileName, removeExcMets, removePrefix, supressWarnings);
%
% Notes
% -----
% A number of consistency checks are performed in order to ensure that the
% model is valid. Take these warnings seriously and modify the model
% structure to solve them.

p=parseRAVENargs(varargin, {'fileName',[]; 'removeExcMets',[]; 'removePrefix',[]; 'supressWarnings',false});
fileName=p.fileName; removeExcMets=p.removeExcMets; removePrefix=p.removePrefix; supressWarnings=p.supressWarnings;
if isempty(fileName)
    [fileName, pathName] = uigetfile({'*.xml;*.sbml'}, 'Please select the model file');
    if fileName == 0
        error('You should select a model file')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);
if isempty(removeExcMets)
    removeExcMets=true;
end

if isempty(removePrefix)
    removePrefix=true;
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

%Set a canonical order of the model fields
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

%RAVEN only supports SBML Level 3 Version 1 with the FBC version 2 package.
%Older SBML levels/versions and non-FBC formats are no longer supported.
if ~(isfield(modelSBML,'SBML_level') && isequal(modelSBML.SBML_level,3) && ...
        isfield(modelSBML,'SBML_version') && isequal(modelSBML.SBML_version,1) && ...
        isfield(modelSBML,'fbc_version') && isequal(modelSBML.fbc_version,2))
    EM=['importModel only supports SBML Level 3 Version 1 with the FBC '...
        'version 2 package. Convert the model to this format first (for '...
        'example with another tool), or use an older RAVEN version to '...
        'import it.'];
    dispEM(EM);
end

%Retrieve compartment names and IDs
compartmentNames=cell(numel(modelSBML.compartment),1);
compartmentIDs=cell(numel(modelSBML.compartment),1);
compartmentOutside=cell(numel(modelSBML.compartment),1);
compartmentMiriams=cell(numel(modelSBML.compartment),1);

if isfield(modelSBML.compartment,'sboTerm') && isscalar(unique([modelSBML.compartment.sboTerm]))
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

%Retrieve info on metabolites
metaboliteNames={};
metaboliteIDs={};
metaboliteCompartments={};
metaboliteUnconstrained=[];
metaboliteFormula={};
metaboliteInChI={};
metaboliteMiriams={};
metaboliteCharges=[];

%Gene information is collected from fbc_geneProduct further below
geneMiriams={};
geneShortNames={};
proteins={};

%In FBC v2 all species are metabolites; genes are stored as fbc_geneProduct
%and are handled separately below.
metSBOs = [];
%Regex of compartment names, later to be used to remove from metabolite
%names if present as suffix.
regexCompNames = ['\s?\[((' strjoin({modelSBML.compartment.name},')|(') '))\]$'];
for i=1:numel(modelSBML.species)
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
                if isscalar(compositionIndexes)
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
        metaboliteFormula{numel(metaboliteFormula)+1,1}='';
        metaboliteMiriams{numel(metaboliteMiriams)+1,1}=[];
    end
    %Get SBO term
    if isfield(modelSBML.species(i),'sboTerm') && ~(modelSBML.species(i).sboTerm==-1)
        metSBOs(end+1,1) = modelSBML.species(i).sboTerm;
    end

    %Remove trailing [compartment] from metabolite name if present
    metaboliteNames{end,1}=regexprep(metaboliteNames{end,1},regexCompNames,'');
    if isfield(modelSBML.species(i),'fbc_charge')
        if ~isempty(modelSBML.species(i).fbc_charge) && modelSBML.species(i).isSetfbc_charge
            metaboliteCharges(numel(metaboliteCharges)+1,1)=double(modelSBML.species(i).fbc_charge);
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

%Add SBO terms to metabolite miriam fields
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

if isfield(modelSBML.reaction,'sboTerm') && isscalar(unique([modelSBML.reaction.sboTerm]))
    %If all the SBO terms are identical, don't add them to rxnMiriams
    modelSBML.reaction = rmfield(modelSBML.reaction,'sboTerm');
end

for i=1:numel(modelSBML.reaction)

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
        [~,fluxBoundIdx] = ismember({lb,ub},parameter.name);
        reactionLB(counter) = parameter.value{fluxBoundIdx(1)};
        reactionUB(counter) = parameter.value{fluxBoundIdx(2)};
    else
        if reactionReversibility(counter)==true
            reactionLB(counter)=-inf;
        else
            reactionLB(counter)=0;
        end
        reactionUB(counter)=inf;
        reactionObjective(counter)=0;
    end

    %Find the associated gene if available. In FBC v2 gene association data
    %is stored in the fbc_geneProductAssociation field.
    if isfield(modelSBML.reaction(i),'fbc_geneProductAssociation')
        if ~isempty(modelSBML.reaction(i).fbc_geneProductAssociation) && ~isempty(modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association)
            grRules{counter}=modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association.fbc_association;
        end
    end

    %Add reaction compartment
    rxnComp='';
    if isfield(modelSBML.reaction(i),'compartment')
        if ~isempty(modelSBML.reaction(i).compartment)
            rxnComp=modelSBML.reaction(i).compartment;
        end
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
        if strfind(modelSBML.reaction(i).annotation,'http://identifiers.org/ec-code')
            eccode=parseAnnotation(modelSBML.reaction(i).annotation,'http://identifiers.org/','/','ec-code');
        elseif strfind(modelSBML.reaction(i).annotation,'https://identifiers.org/ec-code')
            eccode=parseAnnotation(modelSBML.reaction(i).annotation,'https://identifiers.org/','/','ec-code');
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
                idx=ismember(reactionIDs,rxn);
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
    %Finding whether the taxonomy is written in the old or the new way
    J=strfind(modelSBML.annotation,'"http://identifiers.org/');
    if any(J)
        I = I(find(I>J,1))-1;
        J = J(find(J<I,1))+24;
        model.annotation.taxonomy=modelSBML.annotation(J:I);
    else
        J=strfind(modelSBML.annotation,'"https://identifiers.org/');
        if any(J)
            I = I(find(I>J,1))-1;
            J = J(find(J<I,1))+25;
            model.annotation.taxonomy=modelSBML.annotation(J:I);
        end
    end
end
if isfield(modelSBML,'notes')
    startString=strfind(modelSBML.notes,'xhtml">');
    endString=strfind(modelSBML.notes,'</body>');
    if any(startString) && any(endString)
        model.annotation.note=modelSBML.notes(startString(1)+7:endString-1);
        model.annotation.note=regexprep(model.annotation.note,'<p.*?>|</p.*?>','');
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

%If any gene associations have been loaded, build the gene list
if ~all(cellfun(@isempty,grRules))
    %If fbc_geneProduct exists, follow the specified gene order, such
    %that matching geneShortNames in function below will work
    if isfield(modelSBML,'fbc_geneProduct')
        genes={modelSBML.fbc_geneProduct.fbc_id};

        %Get gene Miriams from fbc_geneProduct
        if isempty(geneMiriams)
            geneMiriams = cell(numel(genes),1);
            if isfield(modelSBML.fbc_geneProduct,'sboTerm') && isscalar(unique([modelSBML.fbc_geneProduct.sboTerm]))
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
        genes=getGenesFromGrRules(grRules);
    end
    model.genes=genes;
    model.grRules=grRules;
    [grRules,rxnGeneMat] = standardizeGrRules(model,true);
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;
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

if removePrefix
    [model, hasChanged]=removeIdentifierPrefix(model);
    dispEM(['The following fields have prefixes removed from all entries. '...
    'If this is undesired, run importModel with removePrefix as false. Example: '...
    'importModel(''filename.xml'',[],false);'],false,hasChanged)
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
    % If all subSystems have single entries, then unnest them
    if all(cellfun(@(x) iscell(x) && isscalar(x), model.subSystems))
        model.subSystems = transpose([model.subSystems{:}]);
    end
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
if strfind(searchString,'http://identifiers.org/')
    startString='http://identifiers.org/';
elseif strfind(searchString,'https://identifiers.org/')
    startString='https://identifiers.org/';
else
    miriamStruct=[];
    return;
end

miriamStruct=[];

searchString=regexprep(searchString,'" />','"/>');
[~,targetString] = regexp(searchString,'<rdf:li rdf:resource=".*?"/>','tokens','match');
targetString=regexprep(targetString,'<rdf:li rdf:resource="|"/>','');
targetString=regexprep(targetString,startString,'');

fwdslash  = contains(targetString,'/');
midString = cell(numel(targetString),1);
midString(fwdslash) = {'/'};
midString(~fwdslash) = {':'};

counter=0;
for i=1:numel(targetString)
    if isempty(regexp(targetString{1,i},'inchi|ec-code', 'once'))
        counter=counter+1;
        miriamStruct.name{counter,1} = regexprep(targetString{1,i},[midString{i} '.+'],'','once');
        miriamStruct.value{counter,1} = regexprep(targetString{1,i},[miriamStruct.name{counter,1} midString{i}],'','once');
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
