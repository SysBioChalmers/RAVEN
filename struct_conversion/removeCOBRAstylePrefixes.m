function model=removeCOBRAstylePrefixes(model,fields,forceRemove)
% removeCOBRAstylePrefixes
%   This function removes COBRA style identifier prefixes, which are:
%       "R_" for model.rxns, model.rxnNames and model.id,
%       "M_" for model.mets and model.metNames,
%       "C_" for model.comps;
%       "G_" for model.genes (and also represented in model.grRules).
%   By default, the prefixes are only removed if all entries in a
%   particular field has the prefix.
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
%
% Usage: model=removeCOBRAstylePrefixes(model,fields,forceRemove)

if nargin<2 || isempty(fields)
    fields = {'rxns','mets','comps','genes','metNames','rxnNames','id'};
end
if nargin<3 || isempty(forceRemove)
    forceRemove = true;
end

modelFields = {'rxns',      'R_';
    'mets',      'M_';
    'comps',     'C_';
    'genes',     'G_';
    'metNames',  'M_';
    'rxnNames',  'R_';
    'id',        'M_'};

toChangeIdx = find(ismember(modelFields(:,1),fields));
for i=1:numel(toChangeIdx)
    currName    = modelFields{toChangeIdx(i),1};
    currPrefix  = modelFields{toChangeIdx(i),2};
    currField   = model.(currName);

    if forceRemove
        hasPrefix = true;
    else
        hasPrefix = all(startsWith(currField,currPrefix));
    end
    if hasPrefix
        currField = regexprep(currField,['^' currPrefix],'');
        if strcmp(currName,'genes')
            model.grRules=regexprep(model.grRules,'^G_','');
            model.grRules=regexprep(model.grRules,'\(G_','(');
            model.grRules=regexprep(model.grRules,' G_',' ');
        end
    end
    model.(currName) = currField;
end
end
