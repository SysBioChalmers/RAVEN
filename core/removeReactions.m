function reducedModel=removeReactions(model,rxnsToRemove,removeUnusedMets,removeUnusedGenes,removeUnusedComps)
% removeReactions
%   Deletes a set of reactions from a model
%
%   model             a model structure
%   rxnsToRemove      either a cell array of reaction IDs, a logical vector
%                     with the same number of elements as reactions in the model,
%                     or a vector of indexes to remove
%   removeUnusedMets  remove metabolites that are no longer in use (opt,
%                     default false)
%   removeUnusedGenes remove genes that are no longer in use (opt, default
%                     false)
%   removeUnusedComps remove compartments that are no longer in use (opt,
%                     default false)
%
%   reducedModel      an updated model structure
%
%   Usage: reducedModel=removeReactions(model,rxnsToRemove,removeUnusedMets,...
%           removeUnusedGenes,removeUnusedComps)
%
%   Simonas Marcisauskas, 2016-11-01 - added support for rxnNotes,
%   rxnReferences and confidenceScores
%

if nargin<3
    removeUnusedMets=false;
end
if nargin<4
    removeUnusedGenes=false;
end
if nargin<5
    removeUnusedComps=false;
end

if ischar(rxnsToRemove)
    rxnsToRemove={rxnsToRemove};
end

reducedModel=model;

if ~isempty(rxnsToRemove) || removeUnusedMets || removeUnusedGenes
    indexesToDelete=getIndexes(model,rxnsToRemove,'rxns');

    %Remove reactions
    if ~isempty(indexesToDelete)
        reducedModel.rxns(indexesToDelete)=[];
        if isfield(reducedModel,'lb')
            reducedModel.lb(indexesToDelete)=[];
        end
        if isfield(reducedModel,'ub')
            reducedModel.ub(indexesToDelete)=[];
        end
        if isfield(reducedModel,'rev')
            reducedModel.rev(indexesToDelete)=[];
        end
        if isfield(reducedModel,'c')
            reducedModel.c(indexesToDelete)=[];
        end
        if isfield(reducedModel,'S')
            reducedModel.S(:,indexesToDelete)=[];
        end
        if isfield(reducedModel,'rxnNames')
            reducedModel.rxnNames(indexesToDelete)=[];
        end
        if isfield(reducedModel,'rxnGeneMat')
            reducedModel.rxnGeneMat(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'grRules')
            reducedModel.grRules(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'subSystems')
            reducedModel.subSystems(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'eccodes')
            reducedModel.eccodes(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'equations')
            reducedModel.equations(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnMiriams')
            reducedModel.rxnMiriams(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnComps')
            reducedModel.rxnComps(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnFrom')
            reducedModel.rxnFrom(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnScores')
            reducedModel.rxnScores(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnNotes')
            reducedModel.rxnNotes(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'rxnReferences')
            reducedModel.rxnReferences(indexesToDelete,:)=[];
        end
        if isfield(reducedModel,'confidenceScores')
            reducedModel.confidenceScores(indexesToDelete,:)=[];
        end
    end

    %Remove unused metabolites
    if removeUnusedMets==true
        if isfield(reducedModel,'S')
            [usedMets, ~]=find(reducedModel.S);
            unUsedMets=true(numel(reducedModel.mets),1);
            unUsedMets(usedMets)=false;
            reducedModel=removeMets(reducedModel,unUsedMets,false,false,false,removeUnusedComps);
        end
    end

    %Remove unused genes
    if removeUnusedGenes==true && isfield(reducedModel,'rxnGeneMat')
        %Find all genes that are not used
        [~, b]=find(reducedModel.rxnGeneMat);
        toKeep=false(numel(reducedModel.genes),1);
        toKeep(b)=true;

        reducedModel.genes=reducedModel.genes(toKeep);
        reducedModel.rxnGeneMat=reducedModel.rxnGeneMat(:,toKeep);

        if isfield(reducedModel,'geneShortNames')
           reducedModel.geneShortNames=reducedModel.geneShortNames(toKeep);
        end

        if isfield(reducedModel,'geneMiriams')
           reducedModel.geneMiriams=reducedModel.geneMiriams(toKeep);
        end

        if isfield(reducedModel,'geneFrom')
           reducedModel.geneFrom=reducedModel.geneFrom(toKeep);
        end

        if isfield(reducedModel,'geneComps')
           reducedModel.geneComps=reducedModel.geneComps(toKeep);
        end
    end
else
    reducedModel=model;
end
end
