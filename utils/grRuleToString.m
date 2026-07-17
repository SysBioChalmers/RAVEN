function rule=grRuleToString(tree)
% grRuleToString  Render a grRule syntax tree back to a string.
%
% Writes the standard RAVEN format: " and " / " or " as operators, with
% parentheses added only where precedence requires them, so that
% parseGrRule(grRuleToString(tree)) round-trips to the same structure.
%
% Parameters
% ----------
% tree : struct
%     a tree from parseGrRule. An empty tree yields ''.
%
% Returns
% -------
% rule : char
%     the rendered grRule.
%
% Examples
% --------
%     grRuleToString(parseGrRule('(G1 and G2) or G3'))   % '(G1 and G2) or G3'
%
% See also: parseGrRule, standardizeGrRules

if isempty(tree)
    rule='';
    return
end
rule=nodeToString(tree);
end


function str=nodeToString(node)
switch node.type
    case 'gene'
        str=node.id;
    case 'and'
        % An OR child inside an AND must be bracketed to keep precedence.
        parts=cellfun(@(c) bracketIf(c,'or'),node.children,'UniformOutput',false);
        str=strjoin(parts,' and ');
    case 'or'
        % An AND child inside an OR needs no brackets for correctness, but
        % RAVEN writes complexes bracketed, matching standardizeGrRules and
        % the joinOR helper in removeLowScoreGenes.
        parts=cellfun(@(c) bracketIf(c,'and'),node.children,'UniformOutput',false);
        str=strjoin(parts,' or ');
    otherwise
        error('RAVEN:badGrRule',['Unexpected grRule node type: ' node.type]);
end
end


function str=bracketIf(node,type)
str=nodeToString(node);
if strcmp(node.type,type)
    str=['(' str ')'];
end
end
