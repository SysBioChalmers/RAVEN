function model=editMiriam(model,type,object,miriamName,miriams,keep)
% editMiriam
%   Change MIRIAM annotation fields, one annotation type at the same time.
%
%   Input:
%       model       model structure
%       type        'met', 'rxn', 'gene' or 'comp' dependent on which
%                   objects the annotations should be assigned to
%       object      either a cell array of IDs, a logical vector with the
%                   same number of elements as the type (see above) in the
%                   model, a vector of indexes, or 'all'
%       miriamName  string specifying the namespace of the identifier, for
%                   instance 'bigg.metabolite'. Should be a valid prefix
%                   from identifiers.org (e.g.
%                   https://registry.identifiers.org/registry/bigg.metabolite)
%       miriam      string or cell array of strings with annotation
%                   identifiers, e.g. '12dgr161'
%       keep        one of the following strings, specifying what should be
%                   done if an object already has an existing MIRIAM
%                   annotations with the same miriamName:
%                   'replace'   discard all existing annotations, all will
%                               be overwritten, even if the new annotation
%                               is an empty field. Should only be used if
%                               you do not want to keep any of the old
%                               annotation with the same miriamName
%                   'fill'      only add annotations to those objects that
%                               did not yet have an annotation with that
%                               miriamName. Otherwise, the existing
%                               annotation is kept, even if it is different
%                               from the suggested new annotation
%                   'add'       keep all existing annotations, and add any
%                               new annotations, after removing duplicates
%                   
%   Ouput:
%       model       model structure with updated MIRIAM annotation field
%   
%   Usage: model=editMiriam(model,type,object,miriamName,miriams,keep)
if ~any(strcmp(keep,{'replace','fill','add'}))
    error('Invalid ''keep'', should be ''replace'',''fill'',''add''.')
end
if ~any(strcmp(type,{'met','gene','rxn','comp'}))
    error('Invalid ''type'', should be ''met'', ''gene'', ''rxn'' or ''comp''.')
end
if islogical(object)
    idx=find(object);
elseif isnumeric(object)
    idx=object;
elseif ischar(object)
    if strcmp(object,'all')
        idx=1:numel(model.([type,'s']));
    else
        idx=getIndexes(model,[type,'s'],object);
        if ~all(idx)
            dispEM('The following objects cannot be found in the model: ',true,object(~idx))
        end
    end
else
    [~,idx]=ismember(object,model.([type,'s']));
    if ~all(idx)
        dispEM('The following objects cannot be found in the model: ',true,object(~idx))
    end
end
if ischar(miriams)
    miriams={miriams};
end
if numel(miriams)==1 && numel(idx)~=1
    miriams=repmat(miriams,numel(idx),1);
elseif numel(miriams)~=numel(idx)
    error('The number of annotations does not match the number of objects.')
end

miriamFieldName=[type,'Miriams'];

[extractedMiriams,extractedMiriamNames]=extractMiriam(model.(miriamFieldName),'all');
extractedMiriams=regexprep(extractedMiriams,'^.+\/','');

[~, midx]=ismember(miriamName,extractedMiriamNames);
if strcmp(keep,'replace')
    extractedMiriams(idx,midx)=miriams;
else
    noMiriam=cellfun(@isempty,extractedMiriams(idx,midx));
    extractedMiriams(idx(noMiriam),midx)=miriams(noMiriam);
    if strcmp(keep,'add') % Skipped when only filling empty entries
        strcat(extractedMiriams(idx(~noMiriam),midx), {'; '}, miriams(~noMiriam));
    end
end

% Make Miriam field again
for i=1:numel(model.([type,'s']))
    miriam.name=cell(1,1);
    miriam.value=cell(1,1);
    for j=1:numel(extractedMiriamNames)
        if ~isempty(extractedMiriams{i,j})
            values=strsplit(extractedMiriams{i,j},'; ');
            miriam.name(end+1:end+numel(values),1)=extractedMiriamNames(j);
            miriam.value(end+1:end+numel(values),1)=values;
        end
    end
    miriam.name(1)=[];
    miriam.value(1)=[];
    model.(miriamFieldName){i}=miriam;
end
end
