function produced=canProduce(model,varargin)
% canProduce  Check which metabolites can be produced from a model.
%
% Checks which metabolites can be produced from a model using the
% specified constraints. This is a less advanced but faster version of
% checkProduction.
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
% produced : logical
%     vector with true if the corresponding metabolite could be produced.
%
% Examples
% --------
%     produced = canProduce(model, mets);
%
% See also
% --------
% checkProduction

p=parseRAVENargs(varargin, {'mets',[]});
mets=p.mets;

if isempty(mets)
    mets=model.mets;
elseif ~islogical(mets) && ~isnumeric(mets)
    mets=convertCharArray(mets);
end

[model, rxns]=addExchangeRxns(model,'out',mets);
produced=haveFlux(model,10^-5,rxns);
end
