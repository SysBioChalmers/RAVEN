function canExchange=canExchange(model,direction,varargin)
% canExchange  Check which metabolites a model can produce or consume.
%
% Adds an exchange reaction for each of the specified metabolites and
% checks which of them can carry a flux. This is the unified replacement
% for canProduce ('produce') and canConsume ('consume').
%
% This is a less advanced but faster version of checkProduction, which
% additionally reports which metabolites must be connected in order to make
% the remaining ones producible.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% direction : char
%     'produce' to check for net production of each metabolite, or
%     'consume' to check for net consumption.
%
% Name-Value Arguments
% --------------------
% mets : cell or logical or double
%     either a cell array of metabolite IDs, a logical vector with the same
%     number of elements as metabolites in the model, or a vector of
%     indexes to check for (default model.mets).
%
% Returns
% -------
% canExchange : logical
%     vector with true if the corresponding metabolite could be produced
%     (direction 'produce') or consumed (direction 'consume').
%
% Examples
% --------
%     produced = canExchange(model, 'produce');
%     consumed = canExchange(model, 'consume', mets);
%
% See also
% --------
% checkProduction, findLeakMetabolite

direction=char(direction);
if ~any(strcmp(direction,{'produce','consume'}))
    error('RAVEN:badInput','direction must be ''produce'' or ''consume''');
end

p=parseRAVENargs(varargin, {'mets',[]});
mets=p.mets;

if isempty(mets)
    mets=model.mets;
elseif ~islogical(mets) && ~isnumeric(mets)
    mets=convertCharArray(mets);
end

if strcmp(direction,'produce')
    exchangeDir='out';
else
    exchangeDir='in';
end

[model, rxns]=addExchangeRxns(model,exchangeDir,mets);
canExchange=haveFlux(model,10^-5,rxns);
end
