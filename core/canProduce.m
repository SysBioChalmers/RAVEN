function produced=canProduce(model,mets)
% canProduce
%   Checks which metabolites that can be produced from a model using the
%   specified constraints. This is a less advanced but faster version of
%   checkProduction.
%
%   model       a model structure
%   mets        either a cell array of metabolite IDs, a logical vector 
%               with the same number of elements as metabolites in the model,
%               or a vector of indexes to check for (opt, default model.mets)
%
%   produced    vector with true if the corresponding metabolite could be
%               produced
%
%   Usage: produced=canProduce(model,mets)

if nargin<2
    mets=model.mets;
end

[model, rxns]=addExchangeRxns(model,'out',mets);
produced=haveFlux(model,10^-5,rxns);
end
