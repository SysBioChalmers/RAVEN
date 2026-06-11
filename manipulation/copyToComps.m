function model=copyToComps(model,toComps,varargin)
% copyToComps  Copy reactions to new compartment(s).
%
% Parameters
% ----------
% model : struct
%     a model structure.
% toComps : cell
%     cell array of compartment ids. If there is no match to model.comps
%     then it is added as a new compartment (see compNames and
%     compOutside).
% rxns : cell or logical or double, optional
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of indexes
%     to copy (default model.rxns).
% deleteOriginal : logical, optional
%     true if the original reactions should be removed, making it move the
%     reactions instead (default false).
% compNames : cell, optional
%     cell array of compartment names. Used if new compartments should be
%     added (default toComps).
% compOutside : cell, optional
%     cell array of the id (as in comps) for the compartment surrounding
%     each of the compartments. Used if new compartments should be added
%     (default all {''}).
%
% Returns
% -------
% model : struct
%     an updated model structure.
%
% Examples
% --------
%     model=copyToComps(model,toComps,rxns,deleteOriginal,compNames,compOutside);
%
% Notes
% -----
% New reactions and metabolites will be named as "id_toComps(i)".

p=parseRAVENargs(varargin, {'rxns',[],[]; 'deleteOriginal',false,@emptyOrLogicalScalar; 'compNames',[],[]; 'compOutside','',[]});
rxns=p.rxns;
if isempty(rxns)
    rxns=model.rxns;
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
end
deleteOriginal=p.deleteOriginal;
compNames=p.compNames;
if isempty(compNames)
    compNames=toComps;
else
    compNames=convertCharArray(compNames);
end
compOutside=p.compOutside;
if ~isempty(compOutside)
    compOutside=convertCharArray(compOutside);
    if length(compOutside) ~= length(compNames)
        error('compOutside and compNames should be of equal size.');
    end
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
    if isfield(modelToAdd,'metFrom')
        modelToAdd = rmfield(modelToAdd,'metFrom');
    end
    if isfield(modelToAdd,'rxnFrom')
        modelToAdd = rmfield(modelToAdd,'rxnFrom');
    end
    if isfield(modelToAdd,'geneFrom')
        modelToAdd = rmfield(modelToAdd,'geneFrom');
    end
  
    %Merge the models
    model=mergeModels({model;modelToAdd},'metNames',[],true);
end

if deleteOriginal==true
    model=removeReactions(model,rxns,true,true,true); %Also delete unused compartments
end

model.id=originalID;
model.name=originalName;
end
