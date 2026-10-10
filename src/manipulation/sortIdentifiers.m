function newModel = sortIdentifiers(model)
% sortIdentifiers  Sort model identifiers alphabetically.
%
% Sort reactions, metabolites, genes and compartments alphabetically by
% their identifier.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Returns
% -------
% newModel : struct
%     an updated model structure with alphabetically sorted identifiers.
%
% Examples
% --------
%     newModel = sortIdentifiers(model);

[~,I]=sort(model.rxns);
newModel=permuteModel(model,I,'rxns');
[~,I]=sort(newModel.mets);
newModel=permuteModel(newModel,I,'mets');
if isfield(newModel,'genes')
    [~,I]=sort(newModel.genes);
    newModel=permuteModel(newModel,I,'genes');
end
[~,I]=sort(newModel.comps);
newModel=permuteModel(newModel,I,'comps');
end