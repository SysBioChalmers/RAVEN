function indexes=getIndexes(model, objects, type, returnLogical)
% getIndexes
%   Retrieves the indexes for a list of reactions or metabolites
%
%   Input:
%   model           a model structure
%   objects         either a cell array of IDs, a logical vector with the
%                   same number of elements as metabolites in the model,
%                   of a vector of indexes
%   type            'rxns', 'mets', or 'genes' depending on what to retrieve
%                   'metnames' queries metabolite names, while 'metcomps'
%                   allows to provide specific metabolites and their
%                   compartments in the format metaboliteName[comp]
%   returnLogical   Sets whether to return a logical array or an array with
%                   the indexes (opt, default false)
%
%   Output:
%   indexes         can be a logical array or a double array depending on
%                   the value of returnLogical
%
% 	Usage: indexes=getIndexes(model, objects, type, returnLogical)

if nargin<4
    returnLogical=false;
end

%If the supplied object is a character array, then convert it to a cell
%array
if ischar(objects)
    objects={objects};
end

indexes=[];

switch type
    case 'rxns'
        searchIn=model.rxns;
    case 'mets'
        searchIn=model.mets;
    case 'genes'
        searchIn=model.genes;
    case 'metnames'
        searchIn=model.metNames;
    case 'metcomps'
        % If provided as metaboliteName[comp], then index search is quite
        % different from general approach, and therefore coded separately.
        mets=regexprep(objects,'(^.+)\[(.+)\]$','$1');
        comps=regexprep(objects,'(^.+)\[(.+)\]$','$2');
        indexes=zeros(1,numel(objects));
        for i=1:numel(objects)
            index=find(strcmp(mets(i),model.metNames));
            index=index(strcmp(comps(i),model.comps(model.metComps(index))));
            if ~isempty(index)
                indexes(i)=index;
            else
                error(['Could not find object ' objects{i} ' in the model']);
            end
        end
        indexes=indexes(:);
        if returnLogical==true
            tempIndexes=false(numel(model.mets),1);
            tempIndexes(indexes)=true;
            indexes=tempIndexes;
        end
        return %None of the remaining function needs to run if metcomps
    otherwise
        if containsOct(objects,{'rxns','mets','metnames','metcomps','genes'})
            error('The second and third parameter provided to run getIndexes are likely switched. Note that the second parameter is the object to find the index for, while the third parameter is the object type (''rxns'', ''mets'', ''metnames'', ''metcomps'' or ''genes'').')
        else
            error('Incorrect value of the ''type'' parameter. Allowed values are ''rxns'', ''mets'', ''metnames'', ''metcomps'' or ''genes''.');
        end
end

if iscell(objects)
    for i=1:numel(objects)
        index=find(strcmp(objects(i),searchIn));
        if strcmpi(type,'metnames')
            indexes{i}=index;
        elseif ~isempty(index)
            indexes(i)=index;
        else
            error(['Could not find object ' objects{i} ' in the model']);
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
if iscell(indexes) && length(indexes)==1
    indexes=indexes{1};
end
end
