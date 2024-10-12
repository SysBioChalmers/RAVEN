function model=addCOBRAstylePrefixes(model,fields,printWarnings)
% addCOBRAstylePrefixes
%   This function adds COBRA style identifier prefixes, which are:
%       "R_" for model.rxns,
%       "M_" for model.mets,
%       "C_" for model.comps;
%       "G_" for model.genes (and also represented in model.grRules).
%   If all entries in a field already have the prefix, then no additional
%   prefix will be added.
%
% Input:
%   model           model whose identifiers should be modified
%   fields          cell array with model field names to which the prefix
%                   should be added, possible values: {'rxns', 'mets', 
%                   'comps', 'genes'} (optional, by default prefixes are
%                   added to all listed model fields will be checked).
%   printWarnings   if warnings should be shown (optional, default true).
%
% Output:
%   model           modified model
%
% Usage: model=addCOBRAstylePrefixes(model,fields,printWarnings)

if nargin<2 || isempty(fields)
    fields = {'rxns','mets','comps','genes'};
end
if nargin<3 || isempty(printWarnings)
    printWarnings = true;
end

modelFields = {'rxns','R_';
               'mets','M_';
               'comps','C_';
               'genes','G_'};

toChangeIdx = find(ismember(modelFields(:,1),fields));
for i=1:numel(toChangeIdx)
    currName    = modelFields{toChangeIdx(i),1};
    currPrefix  = modelFields{toChangeIdx(i),2};
    if isfield(model,currName)
        currField   = model.(currName);
    else
        continue;
    end
    if ~all(startsWith(currField,currPrefix))
        currField = strcat(currPrefix, currField);
        if strcmp(currName,'genes')
                model.grRules = regexprep(model.grRules, '(\<[0-9_a-zA-Z])', 'G_$1');
                model.grRules = regexprep(model.grRules, ' G_or ', ' or ');
                model.grRules = regexprep(model.grRules, ' G_and ', ' and ');
        end
        model.(currName) = currField;
    elseif printWarnings
        warning(['All identifiers in "model.' currName '" already start ' ...
            'with "' currPrefix '", no additional prefix added.'])
    end
end
end
