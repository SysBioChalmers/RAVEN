function newIds=generateNewIds(model,type,prefix,quantity,numLength)
% generateNewIds
%   Generates a list of new metabolite or reaction ids, sequentially
%   numbered with a defined prefix. The model is queried for the highest
%   existing number of that type of id.
%
%   model       model structure
%   type        string specifying type of id, 'rxns' or 'mets'
%   prefix      string specifying prefix to be used in all ids. E.g. 's_'
%               or 'r_'.
%   quantity    number of new ids that should be generated (optional, default 1)
%   numLength   length of numerical part of id. E.g. 4 gives ids like
%               r_0001 and 6 gives ids like r_000001. If the prefix is
%               already used in the model, then the model-defined length
%               will be used instead. (optional, default 4)
%
% Usage: newIds=generateNewIds(model,type,prefix,quantity,numLength)
%   
type=char(type);
prefix=char(prefix);

if type=='rxns'
    existingIds=model.rxns;
elseif type=='mets'
    existingIds=model.mets;
else
    error('type should be either ''rxns'' or ''mets''.')
end
if nargin<5
    numLength=4;
end

% Subset only existingIds that have the prefix
existingIds=existingIds(~cellfun(@isempty,regexp(existingIds,['^' prefix])));

if ~isempty(existingIds)
    existingIds=regexprep(existingIds,['^' prefix],'');
    existingIds=sort(existingIds);
    lastId=existingIds{end};
    numLength=length(lastId);
    lastId=str2double(lastId);
else
    lastId=0;
    fprintf(['No ' type ' ids with prefix "' prefix ...
        '" currently exist in the model. The first new id will be "' ...
        [prefix,num2str(1,['%0' num2str(numLength) 'd'])] '"\n'],'%s')
    lastId=0;
end

newIds=cell(quantity,1);

for k=1:quantity
    newIds{k}=[prefix,num2str(k+lastId,['%0' num2str(numLength) 'd'])];
end
end