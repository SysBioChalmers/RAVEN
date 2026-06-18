function [newModel,remGenes] = removeLowScoreGenes(model,geneScores,varargin)
% removeLowScoreGenes  Remove low-scoring genes from model.
%
% This function removes genes from a model based on their scores, a step used
% by the tINIT package. The function recognizes and differentiates between
% isozymes and subunits of an enzyme complex. Genes are removed from each
% grRule, subject to the following conditions:
%
%     1) At least one gene must remain associated with the reaction.
%     2) Genes involved in a complex (joined by ANDs) are not removed.
%
% Parameters
% ----------
% model : struct
%     model structure from which genes are to be removed.
% geneScores : double
%     a vector of scores associated with the model genes. Genes with a
%     positive score will remain in the model, whereas genes with a negative
%     score will try to be removed.
%
%     If all genes associated with a reaction have a negative score, then the
%     least-negative gene will remain; if there is a tie, one will be selected
%     at random.
%
%     If a negative-scoring gene is a subunit in a complex, it will not be
%     removed; however, the entire complex may be removed. See the following
%     example cases:
%
%     - Original: G1 or (G2 and G3 and G4); Negative: G1, G2; New: G2 and G3
%       and G4.
%     - Original: G1 or (G2 and G3) or (G4 and G5); Negative: G1, G2; New: G4
%       and G5 [using default complexScoring].
%     - Original: (G1 and (G2 or G3) and G4); Negative: G2; New: (G1 and G3
%       and G4).
%
% Name-Value Arguments
% --------------------
% isozymeScoring : char
%     method used to combine the scores of multiple isozymes; 'min', 'max',
%     'median', or 'average' (default 'max').
% complexScoring : char
%     method used to combine the scores of enzyme complex subunits: 'min',
%     'max', 'median', or 'average' (default 'min').
%
% Returns
% -------
% newModel : struct
%     model with updated gene, grRules, and rxnGeneMat fields after attempting
%     to remove negative-score genes.
% remGenes : cell
%     a list of negative-score genes that were fully removed from the model.
%     Negative-score genes that were removed from some grRules but not all
%     will not be included in this list.


p=parseRAVENargs(varargin, {'isozymeScoring',[]; 'complexScoring','min'});
isozymeScoring=p.isozymeScoring;
complexScoring=p.complexScoring;
if isempty(isozymeScoring)
    isozymeScoring = 'max';
end

if ~isequal(size(model.genes),size(geneScores))
    error('The dimensions of model genes and geneScores do not match.');
end

% convert logical operators to symbols
grRules = model.grRules;
if any(contains(grRules,{'&','|'}))
    symbolicRules = true;
else
    symbolicRules = false;
end
grRules = regexprep(grRules,' and ',' & ');
grRules = regexprep(grRules,' or ',' | ');

% get unique list of grRules
[uRules,~,rule_ind] = unique(grRules);

% iterate through each rule
newRules = uRules;  %initialize newRules
for i = 1:numel(uRules)
    if isempty(uRules{i}) || ~contains(uRules{i},'|')
        continue
    elseif contains(uRules{i},'&')
        newRules{i} = processComplexRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    else
        newRules{i} = processSimpleRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    end
end

% re-map unique rules to model
newModel = model;
newModel.grRules = newRules(rule_ind);

% restore original logical operator formatting if it was changed
if ~symbolicRules
    newModel.grRules = regexprep(newModel.grRules, ' \| ', ' or ');
    newModel.grRules = regexprep(newModel.grRules, ' & ', ' and ');
end

% regenerate "genes" and "rxnGeneMat" model fields
[genes,rxnGeneMat] = getGenesFromGrRules(newModel.grRules);
newModel.genes = genes;
newModel.rxnGeneMat = rxnGeneMat;

% update other gene-related fields
remInd = ~ismember(model.genes,newModel.genes);
remGenes = model.genes(remInd);

if isfield(newModel,'geneShortNames')
    newModel.geneShortNames(remInd) = [];
end
if isfield(newModel,'proteins')
    newModel.proteins(remInd) = [];
end
if isfield(newModel,'geneMiriams')
    newModel.geneMiriams(remInd) = [];
end
if isfield(newModel,'geneFrom')
    newModel.geneFrom(remInd) = [];
end
if isfield(newModel,'geneComps')
    newModel.geneComps(remInd) = [];
end


end



function [updatedRule,rScore] = processSimpleRule(rule,genes,gScores,isozymeScoring,complexScoring)
% Either score or modify a reaction gene rule containig only ANDs or ORs.
%
% If the rule contains an enzyme complex (all ANDs), the complex will be
% scored based on the score of its subunits. Subunits without a score (NaN)
% will be excluded from the score calculation.
%
% If the rule contains only isozymes (all ORs), the negative-score genes
% will be removed from the rule. Isozymes without a score (NaN) will not be
% removed from the rule. The resuling rule will then be scored.


% get IDs and indices of genes involved in rule
ruleGenes = unique(regexp(rule,'[^&|\(\) ]+','match'));
[~,geneInd] = ismember(ruleGenes,genes);

% rules with one or no genes remain unchanged
if numel(ruleGenes) < 2
    rScore = gScores(geneInd);
    updatedRule = rule;
    return
end

if ~contains(rule,'&')  % rule contains isozymes
    
    scoreMethod = isozymeScoring;
    negInd = gScores(geneInd) < 0;  % NaNs will return false here
    if all(negInd)
        % deterministic least-negative: sort descending; stable sort breaks
        % ties by original gene order (left-to-right in the rule string)
        [~,sortOrder] = sort(gScores(geneInd),'descend');
        maxInd = sortOrder(1);
        updatedRule = ruleGenes{maxInd};
    elseif sum(~negInd) == 1
        updatedRule = ruleGenes{~negInd};
    else
        updatedRule = strjoin(ruleGenes(~negInd),' | ');
        if startsWith(rule,'(')
            updatedRule = ['(',updatedRule,')'];
        end
    end
    
    % update ruleGenes and their indices
    ruleGenes = unique(regexp(updatedRule,'[^&|\(\) ]+','match'));
    [~,geneInd] = ismember(ruleGenes,genes);
    
elseif ~contains(rule,'|')  % rule contains enzyme complex
    scoreMethod = complexScoring;
    updatedRule = rule;
else
    error('This function cannot handle rules with both "OR" and "AND" expressions.');
end

% score rule
switch lower(scoreMethod)
    case 'min'
        rScore = min(gScores(geneInd),[],'omitnan');
    case 'max'
        rScore = max(gScores(geneInd),[],'omitnan');
    case 'median'
        rScore = median(gScores(geneInd),'omitnan');
    case 'average'
        rScore = mean(gScores(geneInd),'omitnan');
end

end



function updatedRule = processComplexRule(rule,genes,gScores,isozymeScoring,complexScoring)
% Update reactions containing both AND and OR expressions via recursive descent.
%
% Walks the boolean parse tree: OR nodes (isozymes) prune negative-score
% branches while keeping at least one; AND nodes (complexes) are scored as
% a unit and may be removed at the enclosing OR level. Nested isozyme groups
% within a complex are simplified in place.
updatedRule = evalGPR(rule, genes, gScores, isozymeScoring, complexScoring);
end



function [rule, score] = evalGPR(expr, genes, gScores, isozymeScoring, complexScoring)
% Recursively evaluate and simplify a GPR expression.
% NaN-scored genes/complexes are not pruned (treated as non-negative).

expr = strtrim(expr);

% Strip redundant outer parentheses
while numel(expr) >= 2 && expr(1) == '('
    closePos = findMatchingClose(expr, 1);
    if closePos == numel(expr)
        expr = strtrim(expr(2:end-1));
    else
        break
    end
end

% OR split: isozyme level
orTerms = splitDepth0(expr, '|');
if numel(orTerms) > 1
    termScores = nan(1, numel(orTerms));
    termRules  = cell(1,  numel(orTerms));
    for k = 1:numel(orTerms)
        [termRules{k}, termScores(k)] = evalGPR(orTerms{k}, genes, gScores, isozymeScoring, complexScoring);
    end
    % NaN-scored terms are kept; negatives pruned when alternatives exist
    negIdx = termScores < 0;  % NaN -> false, so NaN terms are retained
    if all(negIdx)
        % All terms negative: keep the least-negative (deterministic)
        [~, ord] = sort(termScores, 'descend');
        keepIdx  = ord(1);
    else
        keepIdx = find(~negIdx);
    end
    termRules  = termRules(keepIdx);
    termScores = termScores(keepIdx);
    rule  = joinOR(termRules);
    score = applyScoring(termScores, isozymeScoring);
    return
end

% AND split: complex level
andTerms = splitDepth0(expr, '&');
if numel(andTerms) > 1
    termScores = nan(1, numel(andTerms));
    termRules  = cell(1,  numel(andTerms));
    for k = 1:numel(andTerms)
        [termRules{k}, termScores(k)] = evalGPR(andTerms{k}, genes, gScores, isozymeScoring, complexScoring);
    end
    score = applyScoring(termScores, complexScoring);
    rule  = joinAND(termRules);
    return
end

% Leaf: single gene
geneIdx = find(strcmp(genes, expr), 1);
if ~isempty(geneIdx)
    score = gScores(geneIdx);
else
    score = NaN;
end
rule = expr;
end



function terms = splitDepth0(expr, sep)
% Split expr on ' sep ' at bracket depth 0.
terms = {};
depth = 0;
start = 1;
n     = numel(expr);
i     = 1;
while i <= n
    c = expr(i);
    if c == '('
        depth = depth + 1;
        i = i + 1;
    elseif c == ')'
        depth = depth - 1;
        i = i + 1;
    elseif c == sep && depth == 0 && i >= 2 && i <= n-1 && expr(i-1) == ' ' && expr(i+1) == ' '
        terms{end+1} = strtrim(expr(start:i-2));  %#ok<AGROW>
        start = i + 2;
        i     = i + 2;
    else
        i = i + 1;
    end
end
terms{end+1} = strtrim(expr(start:end));
end



function pos = findMatchingClose(expr, openPos)
% Return position of ')' matching '(' at openPos.
depth = 0;
for i = openPos:numel(expr)
    if     expr(i) == '('; depth = depth + 1;
    elseif expr(i) == ')'; depth = depth - 1;
        if depth == 0; pos = i; return; end
    end
end
pos = numel(expr);
end



function rule = joinOR(rules)
% Join rules with ' | '; wrap AND-containing subterms in parens.
n = numel(rules);
if n > 1
    for k = 1:n
        if contains(rules{k}, '&')
            rules{k} = ['(' rules{k} ')'];
        end
    end
end
rule = strjoin(rules, ' | ');
end



function rule = joinAND(rules)
% Join rules with ' & '; wrap OR-containing subterms in parens.
n = numel(rules);
if n > 1
    for k = 1:n
        if contains(rules{k}, '|')
            rules{k} = ['(' rules{k} ')'];
        end
    end
end
rule = strjoin(rules, ' & ');
end



function s = applyScoring(scores, method)
switch lower(method)
    case 'min';     s = min(scores, [], 'omitnan');
    case 'max';     s = max(scores, [], 'omitnan');
    case 'median';  s = median(scores, 'omitnan');
    case 'average'; s = mean(scores, 'omitnan');
    otherwise;      s = max(scores, [], 'omitnan');
end
end



