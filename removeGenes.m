function [reducedModel,notDeleted]=removeGenes(model,genesToRemove, removeUnusedMets, removeRxnsWithComplexes)
% removeGenes
%   Deletes a set of genes and all reactions corresponding to them from a model
%
%   model                     a model structure
%   genesToRemove             either a cell array of gene IDs, a logical vector 
%                             with the same number of elements as genes in the model,
%                             or a vector of indexes to remove
%   removeUnusedMets          remove metabolites that are no longer in use (opt, default
%                             false)
%   removeRxnsWithComplexes   remove reactions that are dependent on a
%                             complex if only part of the complex is in genesToRemove 
%                             (opt, default false)
%
%   reducedModel              an updated model structure
%   notDeleted                a cell array with the genes that couldn't be
%                             deleted. This is an empty cell if removeRxnsWithComplexes 
%                             is true
%
%   Usage: [reducedModel,notDeleted]=removeGenes(model,genesToRemove,
%           removeUnusedMets, removeRxnsWithComplexes)
%
%   Rasmus Agren, 2010-12-16
%

if nargin<3
    removeUnusedMets=false;
end

if nargin<4
    removeRxnsWithComplexes=false;
end

notDeleted={};

reducedModel=model;

if ~isempty(genesToRemove)
    indexesToDelete=getIndexes(reducedModel,genesToRemove,'genes');

    %Find all reactions for these genes
    if ~isempty(indexesToDelete) && isfield(reducedModel,'rxnGeneMat')
        [a crap crap]=find(reducedModel.rxnGeneMat(:,indexesToDelete));
        a=unique(a);
        if removeRxnsWithComplexes==true
            %Delete those reactions even if they use complexes
            reducedModel=removeRxns(reducedModel,a,removeUnusedMets,true);
        else
            %First check that all part of the complex should be removed
            rxnsToRemove=false(numel(a),1);
            
            %Loop through the reactions that could possibly be removed
            for i=1:numel(a)
                [crap geneIndexes crap]=find(reducedModel.rxnGeneMat(a(i),:));
                
                %Check to see if all should be removed
                if numel(geneIndexes)>1
                    if all(ismember(geneIndexes,indexesToDelete))
                        rxnsToRemove(i)=true;
                    else
                        %Don't remove the reaction
                    end
                else
                   %Only one gene, remove the reaction
                   rxnsToRemove(i)=true;
                end
            end
            
            %Get the genes that will be deleted (any involved in the
            %reactions to be deleted
            [crap geneIndexes crap]=find(reducedModel.rxnGeneMat(a(rxnsToRemove),:));
            if ~isempty(setdiff(indexesToDelete,geneIndexes))
                notDeleted=reducedModel.genes(setdiff(indexesToDelete,geneIndexes));
            end
            
            %Delete the reactions
            reducedModel=removeRxns(reducedModel,a(rxnsToRemove),removeUnusedMets,true);
        end
    end
else
    reducedModel=model;
end
end
