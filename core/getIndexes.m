function indexes=getIndexes(model,objects, type, returnLogical)
% getIndexes
%   Retrieves the indexes for a list of reactions or metabolites
%
%   model           a model structure
%   objects         either a cell array of IDs, a logical vector with the 
%                   same number of elements as metabolites in the model,
%                   of a vector of indexes
%   type            'rxns', 'mets', or 'genes' depending on what to retrieve
%   returnLogical   Sets whether to return a logical array or an array with
%                   the indexes (opt, default false)
%
%   indexes         can be a logical array or a double array depending on
%                   the value of returnLogical
%
% 	Usage: indexes=getIndexes(model,objects, type, returnLogical)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<4
    returnLogical=false;
end

%If the supplied object is a character array, then convert it to a cell
%array
if ischar(objects)
    objects={objects};
end

indexes=[];

if strcmpi(type,'rxns')
    searchIn=model.rxns;
else
    if strcmpi(type,'mets')
        searchIn=model.mets;
    else
        if strcmpi(type,'genes')
            searchIn=model.genes;
        else
            dispEM('Incorrect value of the "type" parameter. Allowed values are "rxns", "mets" or "genes"'); 
        end
    end
end

if iscell(objects)
    for i=1:numel(objects)
        index=find(strcmp(objects(i),searchIn),1);
        if ~isempty(index)
            indexes(i)=index;
        else
            dispEM(['Could not find object ' objects{i} ' in the model']);
        end
    end
else
    %Now it's either a logical (or 0/1) array or an array with indexes.
    %We want it to be an array with indexes.
    if all(objects)
        %This gets weird if it's all 1
        indexes=objects;
    else
        indexes=find(objects);
    end
end

if returnLogical==true
    tempIndexes=false(numel(searchIn),1);
    tempIndexes(indexes)=true;
    indexes=tempIndexes;
end

indexes=indexes(:);
