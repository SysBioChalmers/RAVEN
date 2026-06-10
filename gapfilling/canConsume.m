function consumed=canConsume(model,mets)
% canConsume  Check which metabolites can be consumed by a model.
%
% Checks which metabolites can be consumed by a model using the specified
% constraints.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% mets : cell or logical or double, optional
%     either a cell array of metabolite IDs, a logical vector with the same
%     number of elements as metabolites in the model, or a vector of
%     indexes to check for (default model.mets).
%
% Returns
% -------
% consumed : logical
%     vector with true if the corresponding metabolite could be produced.
%
% Examples
% --------
%     consumed = canConsume(model, mets);

if nargin<2
    mets=model.mets;
elseif ~islogical(mets) && ~isnumeric(mets)
    mets=convertCharArray(mets);
end

[model, rxns]=addExchangeRxns(model,'in',mets);
consumed=haveFlux(model,10^-5,rxns);
end
