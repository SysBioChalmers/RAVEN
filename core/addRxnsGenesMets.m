function model=addRxnsGenesMets(model,sourcemodel,rxns,add_gene,rxnNote,confidence)
% addRxnsGenesMets
%   Copies reactions from a source model to a new model, including
%   (new) metabolites and genes
%
%   model           draft model where reactions should be copied to
%   sourcemodel     model where reactions and metabolites are sourced from
%   rxns            cell array with reaction IDs (from source model)
%   add_gene        logical whether genes should be copied or not (opt,
%                   default false)
%   rxnNote         string explaining why reactions were copied to model,
%                   is included as newModel.rxnNotes (opt, default
%                   'Added via addRxnsAndMets()')
%   confidence      string, specifying confidence score for all reactions.
%                   Following doi:10.1038/nprot.2009.203 (opt, default 0)
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
%
%   newModel         an updated model structure
%
% 	This function only works if the draft model and source model follow
%	the same metabolite and compartment naming convention. Metabolites are
%	only matched by metaboliteName[compartment]. Useful if one wants to copy
%	additional reactions from source to draft after getModelFromHomology was
%	used involving the same models.
%
%   Usage: newModel=addRxnsGenesMets(model,sourcemodel,rxns,add_gene,rxnNote,confidence)
%
%   Eduard Kerkhoven, 2017-10-20
%

if nargin<6
    confidence=0;
end
if nargin<5
    rxnNote='Added via addRxnsGenesMets()';
end
if nargin<4
    add_gene=false;
end

%% Obtain indexes of reactions in source model
notNewRxn=rxns(ismember(rxns,model.rxns));
rxns=rxns(~ismember(rxns,model.rxns));
if isempty(rxns)
    throw(MException('','All reactions are already in the model.'));
elseif ~isempty(notNewRxn)
    fprintf('\n The following reactions were already present in the model and will not be added:')
    fprintf(strjoin(notNewRxn,'\n'))
end

rxnIdx=find(ismember(sourcemodel.rxns,rxns)); % Get rxnIDs

%% Add new metabolites
metIdx=find(any(sourcemodel.S(:,rxnIdx),2)); % Get metabolite IDs
% Many of the metabolites in are already in the draft model, so only add the new metabolites

% Match by metNames[metComps]. First make these structures for each model.
model.metCompsN = cellstr(num2str(model.metComps));
map = containers.Map(cellstr(num2str(transpose([1:length(model.comps)]))),model.comps);
model.metCompsN = map.values(model.metCompsN);
model.metCompsN = strcat(model.metNames,'[',model.metCompsN,']');

sourcemodel.metCompsN = cellstr(num2str(sourcemodel.metComps));
map = containers.Map(cellstr(num2str(transpose([1:length(sourcemodel.comps)]))),sourcemodel.comps);
sourcemodel.metCompsN = map.values(sourcemodel.metCompsN);
sourcemodel.metCompsN = strcat(sourcemodel.metNames,'[',sourcemodel.metCompsN,']');

newMetCompsN=sourcemodel.metCompsN(metIdx);
notNewMet=newMetCompsN(ismember(newMetCompsN,model.metCompsN));

if ~isempty(notNewMet)
    fprintf('\n\nThe following metabolites were already present in the model and will not be added:\n')
    fprintf(strjoin(notNewMet,'\n'))
end

metIdx=metIdx(~ismember(sourcemodel.metCompsN(metIdx),model.metCompsN));



if ~isempty(metIdx)
    fprintf('\n\nThe following metabolites will be added to the model:\n')
    fprintf(strjoin(sourcemodel.metCompsN(metIdx),'\n'))    
       
    if isfield(sourcemodel,'mets')
        metsToAdd.mets=sourcemodel.mets(metIdx);
    end
    if isfield(sourcemodel,'metNames')
        metsToAdd.metNames=sourcemodel.metNames(metIdx);
    end
    if isfield(sourcemodel,'metFormulas')
        metsToAdd.metFormulas=sourcemodel.metFormulas(metIdx);
    end
    if isfield(sourcemodel,'metCharge')
        metsToAdd.metCharge=sourcemodel.metCharge(metIdx);
    end
    if isfield(sourcemodel,'metMiriams')
        metsToAdd.metMiriams=sourcemodel.metMiriams(metIdx);
    end
    if isfield(sourcemodel,'metFormulas')
        metsToAdd.metFormulas=sourcemodel.metFormulas(metIdx);
    end
    
    metsToAdd.compartments=strtrim(cellstr(num2str(sourcemodel.metComps(metIdx)))); % Convert from compartment string to comparment number
    [~,idx]=ismember(metsToAdd.compartments,strsplit(num2str(1:length(sourcemodel.comps)))); % Match compartment number to compartment abbreviation
    metsToAdd.compartments=sourcemodel.comps(idx); % Fill in compartment abbreviations

    model=addMets(model,metsToAdd);
end
fprintf('\n\nNumber of metabolites added to the model:\n')
fprintf(num2str(numel(metIdx)))

%% Add new genes
if add_gene
    rxnToAdd.grRules=sourcemodel.grRules(rxnIdx); % Get the relevant grRules
    geneList=strjoin(rxnToAdd.grRules);
    geneList=regexp(geneList,' |)|(|and|or','split');% Remove all grRule punctuation
    geneList=geneList(~cellfun(@isempty,geneList)); % Remove spaces and empty genes
    genesToAdd.genes=setdiff(unique(geneList),model.genes); % Only keep new genes
    if ~isempty(genesToAdd.genes)
        genesToAdd.geneComps=zeros(1,numel(genesToAdd.genes));
        genesToAdd.geneComps(:)=sourcemodel.geneComps(1); % Assume all genes are in same compartment
        model=addGenes(model,genesToAdd);
        fprintf('\n\nNumber of genes added to the model:\n')
        fprintf(num2str(numel(genesToAdd.genes)))
    else
        fprintf('\n\nNo genes added to the model, because no genes were annotated or all genes were already present.')
    end
end

%% Add new reactions
rxnToAdd.equations=constructEquations(sourcemodel,rxnIdx);
rxnToAdd.rxnNames=sourcemodel.rxnNames(rxnIdx);
rxnToAdd.rxns=sourcemodel.rxns(rxnIdx);
rxnToAdd.lb=sourcemodel.lb(rxnIdx);
rxnToAdd.ub=sourcemodel.ub(rxnIdx);
rxnToAdd.rxnNotes=cell(1,numel(rxnToAdd.rxns));
rxnToAdd.rxnNotes(:)={rxnNote};
rxnToAdd.confidenceScores=cell(1,numel(rxnToAdd.rxns));
rxnToAdd.confidenceScores(:)={confidence};
if isfield(sourcemodel,'subSystems')
	rxnToAdd.subSystems=sourcemodel.subSystems(rxnIdx);
end
if isfield(sourcemodel,'eccodes')
	rxnToAdd.eccodes=sourcemodel.eccodes(rxnIdx);
end
model=addRxns(model,rxnToAdd,3,'',false);

fprintf('\n\nNumber of reactions added to the model:\n')
fprintf([num2str(numel(rxnIdx)),'\n'])

end