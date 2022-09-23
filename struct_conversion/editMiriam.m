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
miriamName=char(miriamName);
miriams=convertCharArray(miriams);

%Check 'keep' input
keep=char(keep);
if ~any(strcmp(keep,{'replace','fill','add'}))
    error('Invalid ''keep'', should be ''replace'',''fill'',''add''.')
end
%Check 'type' input
type=char(type);
if ~any(strcmp(type,{'met','gene','rxn','comp'}))
    error('Invalid ''type'', should be ''met'', ''gene'', ''rxn'' or ''comp''.')
end
%Check 'object' input
if islogical(object)
    idxInModel=find(object);
elseif isnumeric(object)
    idxInModel=object;
else
    object=convertCharArray(object);
    if numel(object)==1 && strcmp(object{1},'all')
        idxInModel=1:numel(model.([type,'s']));
    else
        idxInModel=getIndexes(model,object,[type,'s']);
        if ~all(idxInModel)
            dispEM('The following objects cannot be found in the model: ',true,object(~idxInModel))
        end
    end
end
%Check 'miriams' input
if numel(miriams)==1 && numel(idxInModel)~=1
    miriams=repmat(miriams,numel(idxInModel),1);
elseif numel(miriams)~=numel(idxInModel)
    error('The number of annotations does not match the number of objects.')
end

miriamFieldName=[type,'Miriams'];

if isfield(model,miriamFieldName)
    [extractedMiriams,extractedMiriamNames]=extractMiriam(model.(miriamFieldName)(idxInModel),'all');
else
    extractedMiriams=cell(numel(idxInModel),1);
    extractedMiriamNames={miriamName};
end

[~, midx]=ismember(miriamName,extractedMiriamNames);
if midx==0
    extractedMiriamNames(end+1)={miriamName};
    extractedMiriams(:,end+1)=miriams;
elseif strcmp(keep,'replace')
    extractedMiriams(:,midx)=miriams;
else
    noMiriam=cellfun(@isempty,extractedMiriams(:,midx));
    extractedMiriams(noMiriam,midx)=miriams(noMiriam);
    if strcmp(keep,'add') % Skipped when only filling empty entries
        existingMiriams=[split(extractedMiriams(~noMiriam,midx),'; '), miriams(~noMiriam)];
        uniqueMiriams=cell(size(existingMiriams,1),1);
        for i=1:size(existingMiriams,1)
            uniqueMiriams{i,1}=strjoin(unique(existingMiriams(i,:)), {'; '});
        end
        extractedMiriams(~noMiriam,midx)=uniqueMiriams;
    end
end

% Make Miriam field again
for i=1:numel(idxInModel)
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
    model.(miriamFieldName){idxInModel(i),1}=miriam;
end
end
