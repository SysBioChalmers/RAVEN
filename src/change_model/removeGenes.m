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
    [grRules,rxnGeneMat,toCheck] = standardizeGrRules(model,true);
    model.grRules    = grRules;
    model.rxnGeneMat = rxnGeneMat;
else
    toCheck    = [];
    rxnGeneMat = model.rxnGeneMat;
end
reducedModel = model;
%Only remove genes that are actually in the model
try
    ischar(genesToRemove{1});
    genesToRemove=genesToRemove(ismember(genesToRemove,model.genes));
end
if ~isempty(genesToRemove)
    indexesToRemove = getIndexes(model,genesToRemove,'genes');
    if ~isempty(indexesToRemove)
        %Make 0 corresponding columns from rxnGeneMat:
        reducedModel.rxnGeneMat(:,indexesToRemove) = 0;        
        genes = model.genes(indexesToRemove);
        %Check if conflicting grRules contain any of the genes to remove
        if ~isempty(toCheck)
            for i=1:numel(toCheck)
                index   = toCheck(i);
                gIdxs   = find(rxnGeneMat(index,:));
                g2check = model.genes(gIdxs);
                if any(ismember(g2check,genes))
                    warning(['Conflicting grRule #' num2str(index) ' (' model.rxns{index} ') contains at least one of the genes to be removed, this grRule will be bypassed in order to avoid logical errors'])
                end
            end
        end
        canCarryFlux = true(size(model.rxns));       
        %Loop through genes and adapt rxns:
        for i = 1:length(genes)
            %Find all reactions for this gene and loop through them:
            geneRxns = find(rxnGeneMat(:,indexesToRemove(i)));
            if ~isempty(geneRxns)
                for j = 1:numel(geneRxns)
                    index  = geneRxns(j);
                    grRule = reducedModel.grRules{index};
                    ruleGenes = reducedModel.genes(logical(rxnGeneMat(index,:)));
                    if ~ismember(index,toCheck) && canCarryFlux(index) && ~isempty(grRule)
                        %Check if rxn can carry flux without this gene:
                        canCarryFlux(index) = canRxnCarryFlux(ruleGenes,grRule,genes{i});
                        %Adapt gene rule & gene matrix:
                        grRule = removeGeneFromRule(grRule,genes{i});
                        reducedModel.grRules{index} = grRule;
                    end
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
%Format grRules and rxnGeneMatrix after all modifications
if standardizeRules
    [grRules,rxnGeneMat] = standardizeGrRules(reducedModel,true);
    reducedModel.grRules = grRules;
    reducedModel.rxnGeneMat = rxnGeneMat;
end
end

function canIt = canRxnCarryFlux(ruleGenes,geneRule,geneToRemove)
%This function converts a gene rule to a logical statement, and then
%asseses if the rule is true (i.e. rxn can still carry flux) or not (cannot
%carry flux).
geneRule = [' ', geneRule, ' '];
for i = 1:length(ruleGenes)
    if strcmp(ruleGenes{i},geneToRemove)
        geneRule = strrep(geneRule,[' ' ruleGenes{i} ' '],' false ');
        geneRule = strrep(geneRule,['(' ruleGenes{i} ' '],'(false ');
        geneRule = strrep(geneRule,[' ' ruleGenes{i} ')'],' false)');
    else
        geneRule = strrep(geneRule,[' ' ruleGenes{i} ' '],' true ');
        geneRule = strrep(geneRule,['(' ruleGenes{i} ' '],'(true ');
        geneRule = strrep(geneRule,[' ' ruleGenes{i} ')'],' true)');
    end
end
geneRule = strtrim(geneRule);
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
end
