function clauses=grRuleToDNF(rule)
% grRuleToDNF  Convert a grRule to disjunctive normal form.
%
% Returns the rule as a list of AND-clauses, applying distributivity where
% needed: "G1 and (G2 or G3)" becomes {{G1,G2},{G1,G3}}. Each clause is one
% isozyme; the genes within it are the subunits of that complex.
%
% Parameters
% ----------
% rule : char or string or struct
%     a grRule, or a tree from parseGrRule.
%
% Returns
% -------
% clauses : cell
%     one cell array of gene IDs per AND-clause. An empty rule yields {}. A
%     rule with no OR anywhere yields a single clause.
%
% Notes
% -----
% Gene order within a clause follows the order the genes appear in the rule.
% Duplicates are not removed: "G1 and G1" stays a two-element clause, since
% collapsing it would silently rewrite the user's rule.
%
% The number of clauses is the product of the OR-arities across a complex,
% so a pathological rule can expand combinatorially. Callers that build one
% object per clause should check numel(clauses) first.
%
% Examples
% --------
%     grRuleToDNF('(G1 and G2) or G3')     % {{'G1','G2'}, {'G3'}}
%     grRuleToDNF('G1 and (G2 or G3)')     % {{'G1','G2'}, {'G1','G3'}}
%
% See also: parseGrRule, isDnfGrRule, expandModel

if ~isstruct(rule)
    rule=parseGrRule(rule);
end
if isempty(rule)
    clauses={};
    return
end
clauses=nodeToDNF(rule);
end


function clauses=nodeToDNF(node)
% Mirrors _node_to_dnf in raven-toolbox manipulation/expand.py.
switch node.type
    case 'gene'
        clauses={{node.id}};
    case 'or'
        % OR: concatenate the disjuncts' clauses.
        clauses={};
        for i=1:numel(node.children)
            clauses=[clauses nodeToDNF(node.children{i})]; %#ok<AGROW>
        end
    case 'and'
        % AND: cross-product of the children's clauses.
        clauses={{}};
        for i=1:numel(node.children)
            childDNF=nodeToDNF(node.children{i});
            newClauses=cell(1,numel(clauses)*numel(childDNF));
            k=0;
            for a=1:numel(clauses)
                for b=1:numel(childDNF)
                    k=k+1;
                    newClauses{k}=[clauses{a} childDNF{b}];
                end
            end
            clauses=newClauses;
        end
    otherwise
        error('RAVEN:badGrRule',['Unexpected grRule node type: ' node.type]);
end
end
