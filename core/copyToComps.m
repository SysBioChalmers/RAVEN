function model=copyToComps(model,toComps,rxns,deleteOriginal,compNames,compOutside)
% copyToComps
%   Copies reactions to new compartment(s)
%
%   model           a model structure
%   toComps         cell array of compartment ids. If there is no match
%                   to model.comps then it is added as a new compartment
%                   (see below for details)
%   rxns            either a cell array of reaction IDs, a logical vector 
%                   with the same number of elements as reactions in the model,
%                   or a vector of indexes to remove (optional, default
%                   model.rxns)
%   deleteOriginal  true if the original reactions should be removed
%                   (making it move the reactions instead) (optional, default
%                   false)
%   compNames       cell array of compartment names. This is used if new
%                   compartments should be added (optional, default toComps)
%   compOutside     cell array of the id (as in comps) for the compartment
%                   surrounding each of the compartments. This is used if
%                   new compartments should be added (optional, default all {''})
%
%   model           an updated model structure
%
%   NOTE: New reactions and metabolites will be named as "id_toComps(i)".
%
% Usage: model=copyToComps(model,toComps,rxns,deleteOriginal,compNames,compOutside)

if nargin<3
    rxns=model.rxns;
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
end
if nargin<4
    deleteOriginal=false;
end
if nargin<5
    compNames=toComps;
else
    compNames=convertCharArray(compNames);
end
if nargin<6
    compOutside=cell(numel(toComps),1);
    compOutside(:)={''};
else
    compOutside=convertCharArray(compOutside);
end

originalID=model.id;
originalName=model.name;

rxns=getIndexes(model,rxns,'rxns');

for i=1:numel(toComps)
    %Check if the compartment exists, otherwise add it
    [I,J]=ismember(toComps(i),model.comps);
    if I==false
        model.comps=[model.comps;toComps(i)];
        model.compNames=[model.compNames;compNames(i)];
        if isfield(model,'compOutside')
            model.compOutside=[model.compOutside;compOutside(i)];
        end
        if isfield(model,'compMiriams')
            model.compMiriams=[model.compMiriams;cell(1,1)];
        end
        J=numel(model.comps);
    end
    %Copy the reactions by making a model structure with only them, then
    %change the localization, and finally merge with the original model
    modelToAdd=model;
    modelToAdd=removeReactions(modelToAdd,setdiff(1:numel(model.rxns),rxns),true,true);
    modelToAdd.rxns=strcat(modelToAdd.rxns,'_',toComps(i));
    modelToAdd.mets=strcat(modelToAdd.mets,'_',toComps(i));
    modelToAdd.comps=modelToAdd.comps(J);
    modelToAdd.compNames=modelToAdd.compNames(J);
    if isfield(modelToAdd,'compOutside')
        modelToAdd.compOutside=modelToAdd.compOutside(J);
    end
    if isfield(modelToAdd,'compMiriams')
        modelToAdd.compMiriams=modelToAdd.compMiriams(J);
    end
    modelToAdd.metComps=ones(numel(modelToAdd.mets),1);
    
    %Merge the models
    model=mergeModels({model;modelToAdd},'metNames');
end

model=rmfield(model,'rxnFrom');
model=rmfield(model,'metFrom');
model=rmfield(model,'geneFrom');

if deleteOriginal==true
    model=removeReactions(model,rxns,true,true,true); %Also delete unused compartments
end

model.id=originalID;
model.name=originalName;
end
