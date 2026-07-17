function reducedModel = removeGenes(model,genesToRemove,varargin)
% removeGenes  Delete a set of genes from a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% genesToRemove : cell or logical or double
%     either a cell array of gene IDs, a logical vector with the same number
%     of elements as genes in the model, or a vector of indexes to remove.
%
% Name-Value Arguments
% --------------------
% removeUnusedMets : logical
%     remove metabolites that are no longer in use (default false).
% removeBlockedRxns : logical
%     remove reactions that get blocked after deleting the genes (default
%     false).
% standardizeRules : logical
%     format gene rules to be compliant with the standard format (default
%     true).
%
% Returns
% -------
% reducedModel : struct
%     an updated model structure.
%
% Examples
% --------
%     reducedModel = removeGenes(model, genesToRemove, removeUnusedMets, ...
%                        removeBlockedRxns, standardizeRules);

p=parseRAVENargs(varargin, {'removeUnusedMets',false; 'removeBlockedRxns',false; 'standardizeRules',true});
removeUnusedMets=p.removeUnusedMets;
removeBlockedRxns=p.removeBlockedRxns;
standardizeRules=p.standardizeRules;
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
if ~(islogical(genesToRemove) || isnumeric(genesToRemove))
    genesToRemove=convertCharArray(genesToRemove);
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
                    if ~ismember(index,toCheck) && canCarryFlux(index) && ~isempty(grRule)
                        %Check if rxn can carry flux without this gene:
                        canCarryFlux(index) = canRxnCarryFlux(grRule,genes{i});
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

function canIt = canRxnCarryFlux(geneRule,geneToRemove)
%Evaluate the rule with geneToRemove absent and every other gene present.
%A complex (AND) needs all of its subunits, isozymes (OR) need only one.
tree = parseGrRule(geneRule);
if isempty(tree)
    canIt = true;
    return
end
canIt = evalWithout(tree,geneToRemove);
end

function tf = evalWithout(node,geneToRemove)
switch node.type
    case 'gene'
        tf = ~strcmp(node.id,geneToRemove);
    case 'and'
        tf = all(cellfun(@(c) evalWithout(c,geneToRemove),node.children));
    case 'or'
        tf = any(cellfun(@(c) evalWithout(c,geneToRemove),node.children));
    otherwise
        error('RAVEN:badGrRule',['Unexpected grRule node type: ' node.type]);
end
end

function geneRule = removeGeneFromRule(geneRule,geneToRemove)
%This function receives a gene rule and it returns it without the chosen
%gene. A complex that loses a subunit is dropped whole; an isozyme is
%dropped from its OR without disturbing the alternatives.
tree = parseGrRule(geneRule);
geneRule = grRuleToString(pruneGene(tree,geneToRemove));
end

function node = pruneGene(node,geneToRemove)
%Remove geneToRemove from the tree, returning [] when nothing is left that
%can catalyse the reaction.
if isempty(node)
    node = [];
    return
end
switch node.type
    case 'gene'
        if strcmp(node.id,geneToRemove)
            node = [];
        end
    case 'and'
        %A complex missing any subunit cannot form at all.
        for k = 1:numel(node.children)
            if isempty(pruneGene(node.children{k},geneToRemove))
                node = [];
                return
            end
        end
    case 'or'
        kept = {};
        for k = 1:numel(node.children)
            child = pruneGene(node.children{k},geneToRemove);
            if ~isempty(child)
                kept{end+1} = child; %#ok<AGROW>
            end
        end
        if isempty(kept)
            node = [];
        elseif isscalar(kept)
            node = kept{1};
        else
            node.children = kept;
        end
    otherwise
        error('RAVEN:badGrRule',['Unexpected grRule node type: ' node.type]);
end
end
