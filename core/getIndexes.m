function indexes=getIndexes(model, objects, type, returnLogical)
% getIndexes
%   Retrieves the indexes for a list of reactions or metabolites
%
%   model           a model structure
%   objects         either a cell array of IDs, a logical vector with the
%                   same number of elements as metabolites in the model,
%                   of a vector of indexes
%   type            'rxns', 'mets', or 'genes' depending on what to retrieve
%                   'metnames' queries metabolite names, while 'metscomps'
%                   allows to provide specific metabolites and their
%                   compartments in the format metaboliteName[comp]
%   returnLogical   Sets whether to return a logical array or an array with
%                   the indexes (opt, default false)
%
%   indexes         can be a logical array or a double array depending on
%                   the value of returnLogical
%
% 	Usage: indexes=getIndexes(model, objects, type, returnLogical)
%
%   Eduard Kerkhoven, 2018-03-03
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
elseif strcmpi(type,'mets')
    searchIn=model.mets;
elseif strcmpi(type,'genes')
    searchIn=model.genes;
elseif strcmpi(type,'metnames')
    searchIn=model.metNames;
elseif strcmpi(type,'metscomps')
    % If provided as metaboliteName[comp], then index search is quite
    % different from general approach, and therefore coded separately.
    if isstr(objects)
        objects={objects};
    end
    mets=regexprep(objects,'(^.+)\[(.+)\]$','$1');
    comps=regexprep(objects,'(^.+)\[(.+)\]$','$2');
    for i=1:numel(objects)
        index=find(strcmp(mets(i),model.metNames));
        index=index(strcmp(comps(i),model.comps(model.metComps(index))));
        if ~isempty(index)
            indexes(i)=index;
        else
            EM=['Could not find object ' objects{i} ' in the model'];
            dispEM(EM);
        end
    end
    indexes=indexes(:);
    return; % If metscomps is queried, remaining codes doesn't need executing
else
    EM='Incorrect value of the "type" parameter. Allowed values are "rxns", "mets" or "genes"';
    dispEM(EM);
end

if iscell(objects)
    for i=1:numel(objects)
        index=find(strcmp(objects(i),searchIn),1);
        if ~isempty(index)
            indexes(i)=index;
        else
            EM=['Could not find object ' objects{i} ' in the model'];
            dispEM(EM);
        end
    end
else
    %Now it's either a logical (or 0/1) array or an array with indexes. We
    %want it to be an array with indexes.
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
end
