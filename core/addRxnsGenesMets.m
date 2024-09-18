function model=addRxnsGenesMets(model,sourceModel,rxns,addGene,rxnNote,confidence)
% addRxnsGenesMets
%   Copies reactions from a source model to a new model, including
%   (new) metabolites and genes
%
%   model           draft model where reactions should be copied to
%   sourceModel     model where reactions and metabolites are sourced from
%   rxns            cell array with reaction IDs (from source model). Can also
%                   be string if only one reaction is added
%   addGene         three options:
%                   false   no genes are annotated to the new reactions
%                   true    grRules ared copied from the sourceModel and
%                           new genes are added when required
%                   string or cell array
%                           new grRules are specified as string or cell
%                           array, and any new genes are added when
%                           required
%                   (optional, default false)
%   rxnNote         cell array with strings explaining why reactions were copied
%                   to the model, to be included as newModel.rxnNotes. Can also
%                   be string if same rxnNotes should be added for each new
%                   reaction, or only one reaction is to be added (optional, default
%                   'Added via addRxnsAndMets()')
%   confidence      integer specifying confidence score for all reactions.
%                   4:  biochemical data: direct evidence from enzymes
%                       assays
%                   3:  genetic data: knockout/-in or overexpression
%                       analysis
%                   2:  physiological data: indirect evidence, e.g.
%                       secretion products or defined medium requirement
%                       sequence data: genome annotation
%                   1:  modeling data: required for functional model,
%                       hypothetical reaction
%                   0:  no evidence
%                   following doi:10.1038/nprot.2009.203 (optional, default 0)
%
%   newModel        an updated model structure
%
% 	This function only works if the draft model and source model follow
%	the same metabolite and compartment naming convention. Metabolites are
%	only matched by metaboliteName[compartment]. Useful if one wants to copy
%	additional reactions from source to draft after getModelFromHomology was
%	used involving the same models.
%
% Usage: newModel=addRxnsGenesMets(model,sourceModel,rxns,addGene,rxnNote,confidence)

if nargin<6
    confidence=0;
end
rxns=convertCharArray(rxns);
if nargin<5
    rxnNote={'Added via addRxnsGenesMets()'};
else
    rxnNote=convertCharArray(rxnNote);
end
if numel(rxnNote)==1 && numel(rxns)>1
    rxnNoteArray=cell(1,numel(rxns));
    rxnNoteArray(:)=rxnNote;
    rxnNote=rxnNoteArray;
end
if nargin<4
    addGene=false;
end

% Obtain indexes of reactions in source model
[notNewRxn,oldRxn]=ismember(rxns,model.rxns);
rxns=rxns(~notNewRxn);
if isempty(rxns)
    error('All reactions are already in the model.');
elseif ~isempty(notNewRxn)
    fprintf('The following reactions were already present in the model and will not be added:\n')
    fprintf([strjoin(model.rxns(oldRxn(find(oldRxn))),'\n') '\n'])
end

[~, rxnIdx]=ismember(rxns,sourceModel.rxns); % Get rxnIDs
if any(rxnIdx==0)
    dispEM('The following reaction IDs could not be found in the source model:',true,...
        rxns(rxnIdx==0));
end

% Add new metabolites
metIdx=find(any(sourceModel.S(:,rxnIdx),2)); % Get metabolite IDs
% Many of the metabolites in are already in the draft model, so only add
% the new metabolites

% Match by metNames[metComps]. First make these structures for each model.
metCompsN =cellstr(num2str(model.metComps));
map=containers.Map(cellstr(num2str(transpose([1:length(model.comps)]))),model.comps);
metCompsN = map.values(metCompsN);
metCompsN = strcat(model.metNames,'[',metCompsN,']');

sourcemetCompsN=cellstr(num2str(sourceModel.metComps));
map=containers.Map(cellstr(num2str(transpose([1:length(sourceModel.comps)]))),sourceModel.comps);
sourcemetCompsN = map.values(sourcemetCompsN);
sourcemetCompsN=strcat(sourceModel.metNames,'[',sourcemetCompsN,']');

newMetCompsN=sourcemetCompsN(metIdx);
notNewMet=newMetCompsN(ismember(newMetCompsN,metCompsN));

if ~isempty(notNewMet)
    fprintf('\nThe following metabolites were already present in the model and will not be added:\n')
    fprintf([strjoin(transpose(notNewMet),'\n') '\n'])
end

metIdx=metIdx(~ismember(sourcemetCompsN(metIdx),metCompsN));

if ~isempty(metIdx)
    fprintf('\nThe following metabolites will be added to the model:\n')
    fprintf([strjoin(transpose(sourcemetCompsN(metIdx)),'\n') '\n'])
    
    if isfield(sourceModel,'mets')
        metsToAdd.mets=sourceModel.mets(metIdx);
    end
    if isfield(sourceModel,'metNames')
        metsToAdd.metNames=sourceModel.metNames(metIdx);
    end
    if isfield(sourceModel,'metFormulas')
        metsToAdd.metFormulas=sourceModel.metFormulas(metIdx);
    end
    if isfield(sourceModel,'metCharges')
        metsToAdd.metCharges=sourceModel.metCharges(metIdx);
    end
    if isfield(sourceModel,'metMiriams')
        metsToAdd.metMiriams=sourceModel.metMiriams(metIdx);
    end
    if isfield(sourceModel,'metNotes')
        metsToAdd.metNotes=sourceModel.metNotes(metIdx);
    end
    if isfield(sourceModel,'inchis')
        metsToAdd.inchis=sourceModel.inchis(metIdx);
    end
    if isfield(sourceModel,'metSmiles')
        metsToAdd.metSmiles=sourceModel.metSmiles(metIdx);
    end
    if isfield(sourceModel,'metDeltaG')
        metsToAdd.metDeltaG=sourceModel.metDeltaG(metIdx);
    end
    
    metsToAdd.compartments=strtrim(cellstr(num2str(sourceModel.metComps(metIdx)))); % Convert from compartment string to compartment number
    [~,idx]=ismember(metsToAdd.compartments,strsplit(num2str(1:length(sourceModel.comps)))); % Match compartment number to compartment abbreviation
    metsToAdd.compartments=sourceModel.comps(idx); % Fill in compartment abbreviations
    
    model=addMets(model,metsToAdd);
end
fprintf('\nNumber of metabolites added to the model:\n')
fprintf([num2str(numel(metIdx)),'\n'])

% Add new genes
if ~islogical(addGene) | addGene ~= false
    if ischar(addGene)
        rxnToAdd.grRules={addGene};
    elseif iscell(addGene)
        rxnToAdd.grRules=addGene(~notNewRxn); 
    else
        rxnToAdd.grRules=sourceModel.grRules(rxnIdx); % Get the relevant grRules
    end
end
% Add new reactions
rxnToAdd.equations=constructEquations(sourceModel,rxnIdx);
rxnToAdd.rxnNames=sourceModel.rxnNames(rxnIdx);
rxnToAdd.rxns=sourceModel.rxns(rxnIdx);
rxnToAdd.lb=sourceModel.lb(rxnIdx);
rxnToAdd.ub=sourceModel.ub(rxnIdx);
rxnToAdd.rxnNotes(:)=rxnNote(~notNewRxn);
rxnToAdd.rxnConfidenceScores=NaN(1,numel(rxnToAdd.rxns));
if ~isnumeric(confidence)
    EM='confidence score must be numeric';
    dispEM(EM, true);
end
rxnToAdd.rxnConfidenceScores(:)=confidence;
if isfield(sourceModel,'rxnDeltaG')
    rxnToAdd.rxnDeltaG=sourceModel.rxnDeltaG(rxnIdx);
end
if isfield(sourceModel,'subSystems')
    rxnToAdd.subSystems=sourceModel.subSystems(rxnIdx);
end
if isfield(sourceModel,'eccodes')
    rxnToAdd.eccodes=sourceModel.eccodes(rxnIdx);
end
model=addRxns(model,rxnToAdd,3,'',false,true);

fprintf('\nNumber of reactions added to the model:\n')
fprintf([num2str(numel(rxnIdx)),'\n'])
end
