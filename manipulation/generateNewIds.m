function newIds=generateNewIds(model,type,prefix,varargin)
% generateNewIds  Generate a list of new metabolite or reaction ids.
%
% The ids are sequentially numbered with a defined prefix. The model is
% queried for the highest existing number of that type of id.
%
% Parameters
% ----------
% model : struct
%     model structure.
% type : char
%     type of id, 'rxns' or 'mets'.
% prefix : char
%     prefix to be used in all ids, e.g. 's_' or 'r_'.
%
% Name-Value Arguments
% --------------------
% quantity : double
%     number of new ids that should be generated (default 1).
% numLength : double
%     length of the numerical part of the id. E.g. 4 gives ids like r_0001
%     and 6 gives ids like r_000001. If the prefix is already used in the
%     model, then the model-defined length will be used instead
%     (default 4).
%
% Returns
% -------
% newIds : cell
%     cell array with the generated ids.
%
% Examples
% --------
%     newIds = generateNewIds(model, type, prefix, quantity, numLength);
type=char(type);
prefix=char(prefix);

if type=='rxns'
    existingIds=model.rxns;
elseif type=='mets'
    existingIds=model.mets;
else
    error('type should be either ''rxns'' or ''mets''.')
end
p=parseRAVENargs(varargin, {'quantity',1; 'numLength',4});
quantity=p.quantity;
numLength=p.numLength;

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
end

if isnan(lastId)
    lastId = 0;
end

newIds=cell(quantity,1);

for k=1:quantity
    newIds{k}=[prefix,num2str(k+lastId,['%0' num2str(numLength) 'd'])];
end
end