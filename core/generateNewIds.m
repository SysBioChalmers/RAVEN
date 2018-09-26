function newIds=generateNewIds(model,type,prefix,quantity)
% generateNewIds
%   Generates a list of new metabolite or reaction ids, sequentially
%   numbered with a defined prefix. The model is queried for the highest
%   existing number of that type of id.
%
%   model       model structure
%   type        string specifying type of id, 'rxns' or 'mets'
%   prefix      string specifying prefix to be used in all ids. E.g. 's_'
%               or 'r_'.
%   quantity    number of new ids that should be generated
%
%   Usage: newIds=generateNewIds(model,type,prefix,quantity)
%   
%   Eduard Kerkhoven, 2018-09-26
%

%TODO: If prefix is not used yet in model, start numbering from 0001.

if type=='rxns'
    existingIds=model.rxns;
elseif type=='mets'
    existingIds=model.mets;
else
    error('type should be either ''rxns'' or ''mets''.')
end

% Subset only existingIds that have the prefix
existingIds=existingIds(~cellfun(@isempty,regexp(existingIds,['^' prefix])));

existingIds=regexprep(existingIds,['^' prefix],'');
existingIds=sort(existingIds);
lastId=existingIds{end};
idLength=length(lastId);
lastId=str2double(lastId);

newIds=cell(quantity,1);

for k=1:quantity
    newIds{k}=[prefix,num2str(k+lastId,['%0' num2str(idLength) 'd'])];
end
end