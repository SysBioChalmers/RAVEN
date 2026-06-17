function issues=findPotentialErrors(model)
% findPotentialErrors  Identify reactions with non-DNF GPR rules.
%
% Detects grRules that are not in disjunctive normal form (DNF) —
% i.e. rules where an OR operator is nested inside an AND context.
% Standard RAVEN format requires DNF: groups of genes that must all be
% present together (enzyme complexes) are written as conjunctions, and
% the reaction is catalysed by at least one complex (disjunction of
% conjunctions).
%
% Parameters
% ----------
% model : struct
%     a model structure with grRules and rxns fields.
%
% Returns
% -------
% issues : struct
%     struct array — one element per problematic reaction — with fields:
%
%     - index  : double — 1-based index into model.rxns.
%     - rxn    : char   — reaction ID.
%     - grRule : char   — the problematic grRule string.
%     - reason : char   — short explanation.
%
% Notes
% -----
% Detection uses a parse-tree walk rather than substring search, so gene
% IDs that contain substrings such as "and" or "or" do not cause false
% positives or negatives.
%
% A rule is non-DNF when the top-level operator is AND and one of its
% operands contains OR, or equivalently when any OR subexpression appears
% inside a bracket group that is joined to another group by AND.
% Use standardizeGrRules to attempt automatic repair.
%
% Examples
% --------
%     issues = findPotentialErrors(model);
%     if ~isempty(issues)
%         fprintf('%d reactions with non-DNF grRules\n', numel(issues));
%     end
%     indexes2check = vertcat(issues.index);

issues=struct('index',{},'rxn',{},'grRule',{},'reason',{});

if ~isfield(model,'grRules') || ~isfield(model,'rxns')
    return
end

for i=1:numel(model.grRules)
    rule=strtrim(model.grRules{i});
    if isempty(rule)
        continue
    end
    % Normalize operator keywords to single-character symbols so that gene
    % IDs containing " and " or " or " as substrings are not misread.
    ruleN=lower(rule);
    ruleN=regexprep(ruleN,' and ',' & ');
    ruleN=regexprep(ruleN,' or ',' | ');

    % Rules with only one operator type are always DNF.
    hasAnd=contains(ruleN,' & ');
    hasOr =contains(ruleN,' | ');
    if ~hasAnd || ~hasOr
        continue
    end

    % Split at top-level | (OR) to get each conjunction term.
    terms=splitTopLevel(ruleN,'|');
    nonDnf=false;
    for j=1:numel(terms)
        term=strtrim(terms{j});
        % Strip a single pair of enclosing parens if present.
        if ~isempty(term) && term(1)=='(' && matchClose(term,1)==numel(term)
            term=term(2:end-1);
        end
        % If this AND-term still contains an OR operator the rule is non-DNF.
        if contains(term,' | ')
            nonDnf=true;
            break
        end
    end

    if nonDnf
        issues(end+1,1)=struct( ...
            'index',  i, ...
            'rxn',    model.rxns{i}, ...
            'grRule', rule, ...
            'reason', 'Non-DNF: OR operator nested inside AND context'); %#ok<AGROW>
    end
end
end


function parts=splitTopLevel(str,sep)
% Split str on sep only at bracket depth zero.
parts={};
depth=0;
start=1;
n=numel(str);
sn=numel(sep);
k=1;
while k<=n
    c=str(k);
    if c=='('
        depth=depth+1;
    elseif c==')'
        depth=depth-1;
    elseif depth==0 && k+sn-1<=n && strcmp(str(k:k+sn-1),sep)
        parts{end+1}=str(start:k-1); %#ok<AGROW>
        start=k+sn;
        k=k+sn-1;
    end
    k=k+1;
end
parts{end+1}=str(start:end);
end


function pos=matchClose(str,openPos)
% Return the position of the closing ')' that matches str(openPos).
depth=0;
pos=-1;
for k=openPos:numel(str)
    if str(k)=='('
        depth=depth+1;
    elseif str(k)==')'
        depth=depth-1;
        if depth==0
            pos=k;
            return
        end
    end
end
end
