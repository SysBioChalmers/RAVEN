function reducedModel=removeReactions(model,rxnsToRemove,varargin)
% removeReactions  Delete a set of reactions from a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxnsToRemove : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of indexes
%     to remove.
%
% Name-Value Arguments
% --------------------
% removeUnusedMets : logical
%     remove metabolites that are no longer in use (default false).
% removeUnusedGenes : logical
%     remove genes that are no longer in use (default false).
% removeUnusedComps : logical
%     remove compartments that are no longer in use (default false).
%
% Returns
% -------
% reducedModel : struct
%     an updated model structure.
%
% Examples
% --------
%     reducedModel = removeReactions(model, rxnsToRemove, removeUnusedMets, ...
%         removeUnusedGenes, removeUnusedComps);

p=parseRAVENargs(varargin, {'removeUnusedMets',false; 'removeUnusedGenes',false; 'removeUnusedComps',false});
removeUnusedMets=p.removeUnusedMets;
removeUnusedGenes=p.removeUnusedGenes;
removeUnusedComps=p.removeUnusedComps;
if ~islogical(rxnsToRemove) && ~isnumeric(rxnsToRemove)
    rxnsToRemove=convertCharArray(rxnsToRemove);
end

reducedModel=model;
reg=ravenModelFields();

if ~isempty(rxnsToRemove) || removeUnusedMets || removeUnusedGenes
    indexesToDelete=getIndexes(model,rxnsToRemove,'rxns');
    
    %Remove reactions
    if ~isempty(indexesToDelete)
        rxnReg=reg(strcmp({reg.type},'rxn'));
        for ri_=1:numel(rxnReg)
            fname=rxnReg(ri_).name;
            if isfield(reducedModel,fname)
                reducedModel.(fname)(indexesToDelete)=[];
            end
        end
        if isfield(reducedModel,'S')
            reducedModel.S(:,indexesToDelete)=[];
        end
        if isfield(reducedModel,'rxnGeneMat')
            reducedModel.rxnGeneMat(indexesToDelete,:)=[];
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
        
        geneReg=reg(strcmp({reg.type},'gene'));
        for ri_=1:numel(geneReg)
            fname=geneReg(ri_).name;
            if isfield(reducedModel,fname)
                reducedModel.(fname)=reducedModel.(fname)(toKeep);
            end
        end
        reducedModel.rxnGeneMat=reducedModel.rxnGeneMat(:,toKeep);
    end
else
    reducedModel=model;
end
end
