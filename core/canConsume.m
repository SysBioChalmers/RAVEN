function consumed=canConsume(model,mets)
% canConsume
%   Checks which metabolites that can be consumed by a model using the
%   specified constraints
%
%   model       a model structure
%   mets        either a cell array of metabolite IDs, a logical vector 
%               with the same number of elements as metabolites in the model,
%               or a vector of indexes to check for (opt, default model.mets)
%
%   consumed    vector with true if the corresponding metabolite could be
%               produced
%
%   Usage: consumed=canConsume(model,mets)

if nargin<2
    mets=model.mets;
end

[model, rxns]=addExchangeRxns(model,'in',mets);
consumed=haveFlux(model,10^-5,rxns);
end
