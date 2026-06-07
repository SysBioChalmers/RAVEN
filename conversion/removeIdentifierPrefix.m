function [model, hasChanged]=removeIdentifierPrefix(model,fields,forceRemove)
% removeIdentifierPrefix
%   This function removes identifier prefixes:
%       "R_" for model.rxns, model.rxnNames and model.id,
%       "M_" for model.mets and model.metNames,
%       "C_" for model.comps;
%       "G_" for model.genes (and also represented in model.grRules).
%   By default, the prefixes are only removed if all entries in a
%   particular field has the prefix. The prefixes might have been present
%   because one or more identifiers do not start with a letter or _, which
%   conflicts with SBML specifications.
%
% Input:
%   model           model whose identifiers should be modified
%   fields          cell array with model field names from which the
%                   identifiers should be removed, possible values: 
%                   'rxns', 'mets', 'comps', 'genes', 'metNames', 
%                   'rxnNames', 'id'. (optional, by default all listed
%                   model fields will be checked).
%   forceRemove     if prefixes should be removed even if not all entries
%                   in a model field have the prefix (optional, default
%                   false)
%
% Output:
%   model           modified model
%   hasChanged      cell array with fields and prefixes that are removed
%
% Usage: model=removeIdentifierPrefix(model,fields,forceRemove)

if nargin<2 || isempty(fields)
    fields = {'rxns','mets','comps','genes','metNames','rxnNames','id'};
end
if nargin<3 || isempty(forceRemove)
    forceRemove = false;
end

modelFields = {'rxns',      'R_';
    'mets',      'M_';
    'comps',     'C_';
    'genes',     'G_';
    'metNames',  'M_';
    'rxnNames',  'R_';
    'id',        'M_'};

toChangeIdx = find(ismember(modelFields(:,1),fields));
hasChanged  = false(numel(modelFields(:,1)),1);
for i=1:numel(toChangeIdx)
    currName    = modelFields{toChangeIdx(i),1};
    currPrefix  = modelFields{toChangeIdx(i),2};
    currField   = model.(currName);

    if forceRemove && any(startsWith(currField,currPrefix))
        hasPrefix = true;
    else
        hasPrefix = all(startsWith(currField,currPrefix));
    end
    if hasPrefix
        currField = regexprep(currField,['^' currPrefix],'');
        hasChanged(toChangeIdx(i)) = true;
        if strcmp(currName,'genes')
            model.grRules=regexprep(model.grRules,'^G_','');
            model.grRules=regexprep(model.grRules,'\(G_','(');
            model.grRules=regexprep(model.grRules,' G_',' ');
        end
    end
    model.(currName) = currField;
end
hasChanged = modelFields(hasChanged,:);
hasChanged = append('model.', hasChanged(:,1), ' (', hasChanged(:,2), ' prefix)');
end
