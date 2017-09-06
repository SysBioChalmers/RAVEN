function model=addsRxnsGenesMets(model,sourcemodel,rxns,add_gene,rxnNote,confidence)
% addRxnsGenesMets
%   Adds reactions to a model, including new metabolites and genes
%
%   model           draft model where reactions should be added
%   sourcemodel     model where reactions and metabolites are sourced from
%   rxns            cell array with reaction IDs (from source model)
%   add_gene        logical whether genes should be added or not (opt,
%                   default false)
%   rxnNote         string explaining why reactions were added to model,
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
%	only matched on their metabolite ID. Useful if one wants to move
%	additional reactions from source to draft after getModelFromHomology was
%	used involving the same models.
%
%   Usage: newModel=addRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
%
%   Simonas Marcisauskas, 2017-08-25
%

if nargin<6
    confidence=0;
end
if nargin<5
    rxnNote='Added via addRxnsAndMets()';
end
if nargin<4
    add_gene=false;
end


%% Obtain indexes of reactions in source model
notNewRxn=rxns(ismember(rxns,model.rxns));
rxns=rxns(~ismember(rxns,model.rxns));
if isempty(rxns)
    throw(MException('','All reactions are already in the model.'));
end

if ~isempty(notNewRxn)
    disp('The following reactions were already present in the model and will not be added:')
    disp(strjoin(notNewRxn,'\n'))
end


rxnIdx=find(ismember(sourcemodel.rxns,rxns)); % Get rxnIDs

%% Add new metabolites
metIdx=find(any(sourcemodel.S(:,rxnIdx),2)); % Get metabolite IDs
% Many of the metabolites in are already in the draft model, so only add the new metabolites
mets=sourcemodel.mets(metIdx);
notNewMet=mets(ismember(mets,model.mets));
if ~isempty(notNewMet)
    disp('The following metabolites were already present in the model and will not be added:')
    disp(strjoin(notNewMet,'\n'))
end

metIdx=metIdx(~ismember(sourcemodel.mets(metIdx),model.mets));

if ~isempty(metIdx)
    metsToAdd.mets=sourcemodel.mets(metIdx);
    metsToAdd.metNames=sourcemodel.metNames(metIdx);
    metsToAdd.metFormulas=sourcemodel.metFormulas(metIdx);

    metsToAdd.compartments=strtrim(cellstr(num2str(sourcemodel.metComps(metIdx)))); % Convert from compartment string to comparment number
    [~,idx]=ismember(metsToAdd.compartments,strsplit(num2str(1:length(sourcemodel.comps)))); % Match compartment number to compartment abbreviation
    metsToAdd.compartments=sourcemodel.comps(idx); % Fill in compartment abbreviations

    model=addMets(model,metsToAdd);
end
disp('Number of metabolites added to the model:')
disp(numel(metIdx))

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
        disp('Number of genes added to the model:')
        disp(numel(genesToAdd.genes))
    else
        disp(['No genes added to the model, because no genes were annotated or all genes were already present'])
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
rxnToAdd.rxnConfidenceScores=cell(1,numel(rxnToAdd.rxns));
rxnToAdd.rxnConfidenceScores(:)={confidence};
model=addRxns(model,rxnToAdd,3,'',false);

disp('Number of reactions added to the model:')
disp(numel(rxnIdx))

end
