function reducedModel = removeGenes(model,genesToRemove,removeUnusedMets,removeBlockedRxns,standardizeRules)
% removeGenes
%   Deletes a set of genes from a model
%
%   model                   a model structure
%   genesToRemove           either a cell array of gene IDs, a logical vector
%                           with the same number of elements as genes in the model,
%                           or a vector of indexes to remove
%   removeUnusedMets        remove metabolites that are no longer in use (opt, default
%                           false)
%   removeBlockedRxns       remove reactions that get blocked after deleting the genes
%                           (opt, default false)
%   standardizeRules        format gene rules to be compliant with standard format
%                           (opt, default true)
%
%   reducedModel            an updated model structure
%
%   Usage: reducedModel = removeGenes(model,genesToRemove,removeUnusedMets,removeBlockedRxns)
%
%   Benjamín J. Sánchez, 2018-03-29
%

if nargin<3
    removeUnusedMets = false;
end

if nargin<4
    removeBlockedRxns = false;
end

if nargin<5
    standardizeRules = true;
end

%Format grRules and rxnGeneMatrix:
if standardizeRules
    model = standardizeGrRules(model);
end

reducedModel = model;

if ~isempty(genesToRemove)
    indexesToRemove = getIndexes(model,genesToRemove,'genes');
    if ~isempty(indexesToRemove)
        %Make 0 corresponding columns from rxnGeneMat:
        reducedModel.rxnGeneMat(:,indexesToRemove) = 0;
        
        genes        = model.genes(indexesToRemove);
        canCarryFlux = true(size(model.rxns));
        
        %Loop through genes and adapt rxns:
        for i = 1:length(genes)
            %Find all reactions for this gene and loop through them:
            isGeneInRxn = ~cellfun(@isempty,strfind(reducedModel.grRules,genes{i}));
            for j = 1:length(reducedModel.grRules)
                if isGeneInRxn(j) && canCarryFlux(j)
                    grRule = reducedModel.grRules{j};
                    
                    %Check if rxn can carry flux without this gene:
                    canCarryFlux(j) = canRxnCarryFlux(reducedModel,grRule,genes{i});
                    
                    %Adapt gene rule & gene matrix:
                    grRule = removeGeneFromRule(grRule,genes{i});
                    reducedModel.grRules{j} = grRule;
                end
            end
        end
        
        %Delete or block the reactions that cannot carry flux:
        if removeBlockedRxns
            rxnsToRemove = reducedModel.rxns(~canCarryFlux);
            reducedModel = removeReactions(reducedModel,rxnsToRemove,removeUnusedMets,true);
        else
            reducedModel = removeReactions(reducedModel,[],removeUnusedMets,true);
            reducedModel.lb(~canCarryFlux) = 0;
            reducedModel.ub(~canCarryFlux) = 0;
        end
    end
end

end

function canIt = canRxnCarryFlux(model,geneRule,geneToRemove)
    %This function converts a gene rule to a logical statement, and then asseses
    % if the rule is true (i.e. rxn can still carry flux) or not (cannot carry flux).
    for i = 1:length(model.genes)
        if strcmp(model.genes{i},geneToRemove)
            geneRule = strrep(geneRule,model.genes{i},'false');
        else
            geneRule = strrep(geneRule,model.genes{i},'true');
        end
    end
    geneRule = strrep(geneRule,'and','&&');
    geneRule = strrep(geneRule,'or','||');
    canIt    = eval(geneRule);
end

function geneRule = removeGeneFromRule(geneRule,geneToRemove)
    %This function receives a standard gene rule and it returns it without the
    %chosen gene.
    geneSets = strsplit(geneRule,' or ');
    hasGene  = ~cellfun(@isempty,strfind(geneSets,geneToRemove));
    geneSets = geneSets(~hasGene);
    geneRule = strjoin(geneSets,' or ');
    if length(geneSets) == 1 && ~isempty(strfind(geneRule,'('))
        geneRule = geneRule(2:end-1);
    end
end




