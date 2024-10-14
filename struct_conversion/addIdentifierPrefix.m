function [model, hasChanged]=addIdentifierPrefix(model,fields)
% addIdentifierPrefix
%   If reaction, metabolite, compartment, gene or model identifiers do not
%   start with a letter or _, which conflicts with SBML specifications,
%   prefixes are added for all identifiers in the respective model field.
%   The prefixes are:
%       "R_" for model.rxns,
%       "M_" for model.mets,
%       "C_" for model.comps;
%       "G_" for model.genes (and also represented in model.grRules)
%
% Input:
%   model           model whose identifiers should be modified
%   fields          cell array with model field names that should be
%                   checked if prefixes should be added, possible values: 
%                   'rxns', 'mets', 'comps', 'genes', 'id'. (optional, by
%                   default all listed model fields will be checked).
%
% Output:
%   model           modified model
%   hasChanged      cell array with fields and prefixes that are added
%
% Usage: [model, hasChanged]=addIdentifierPrefix(model,fields)

if nargin<2 || isempty(fields)
    fields = {'rxns','mets','comps','genes','id'};
end

modelFields = {'rxns','R_';
               'mets','M_';
               'comps','C_';
               'genes','G_';
               'id','M_'};

toChangeIdx = find(ismember(modelFields(:,1),fields));
hasChanged  = false(numel(modelFields(:,1)),1);
for i=1:numel(toChangeIdx)
    currName    = modelFields{toChangeIdx(i),1};
    currPrefix  = modelFields{toChangeIdx(i),2};
    if isfield(model,currName)
        currField   = model.(currName);
    else
        continue;
    end
    if ~all(startsWith(currField,regexpPattern('^[a-zA-Z_]')))
        currField = strcat(currPrefix, currField);
        hasChanged(toChangeIdx(i)) = true;

        if strcmp(currName,'genes')
                model.grRules = regexprep(model.grRules, '(\<[0-9_a-zA-Z])', 'G_$1');
                model.grRules = regexprep(model.grRules, ' G_or ', ' or ');
                model.grRules = regexprep(model.grRules, ' G_and ', ' and ');
        end
        model.(currName) = currField;
    end
end

hasChanged = modelFields(hasChanged,:);
hasChanged = append('model.', hasChanged(:,1), ' (', hasChanged(:,2), ' prefix)');
end
