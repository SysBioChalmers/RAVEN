function [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,reactionType)
% getExchangeRxns
%   Retrieves the exchange reactions from a model. Exchange reactions are
%   identified by having either no substrates or products.
%
% Input:
%   model               a model structure
%   reactionType        which exchange reactions should be returned
%                       'all'     all reactions, irrespective of reaction
%                                 bounds
%                       'uptake'  reactions with bounds that imply that
%                                 only uptake are allowed. Reaction
%                                 direction, upper and lower bounds are
%                                 all considered
%                       'excrete' reactions with bounds that imply that
%                                 only excretion are allowed. Reaction
%                                 direction, upper and lower bounds are
%                                 all considered
%                       'reverse' reactions with non-zero upper and lower
%                                 bounds that imply that both uptake and
%                                 excretion are allowed
%                       'blocked' reactions that have zero upper and lower
%                                 bounds, not allowing any flux
%                       'in'      reactions where the boundary metabolite
%                                 is the substrate of the reaction, a
%                                 positive flux value would imply uptake,
%                                 but reaction bounds are not considered
%                       'out'     reactions where the boundary metabolite
%                                 is the substrate of the reaction, a
%                                 positive flux value would imply uptake,
%                                 but reaction bounds are not considered.
%
% Output:
%   exchangeRxns        cell array with the IDs of the exchange reactions
%   exchangeRxnsIndexes vector with the indexes of the exchange reactions
%
% Note:
%   The union of 'in' and 'out' equals 'all'. Also, the union of 'uptake',
%   'excrete', 'reverse' and 'blocked' equals all.
%
% Usage: [exchangeRxns,exchangeRxnsIndexes]=getExchangeRxns(model,reactionType)

if nargin<2
    reactionType='all';
else
    reactionType=char(reactionType);
end

% Find exchange reactions
if isfield(model, 'unconstrained')
    [~, I]=find(model.S(model.unconstrained~=0,:)>0);
    hasNoProd(I)=true;
    [~, I]=find(model.S(model.unconstrained~=0,:)<0);
    hasNoSubs(I)=true;
else
    hasNoProd = transpose(find(sum(model.S>0)==0));
    hasNoSubs = transpose(find(sum(model.S<0)==0));
end
allExch   = [hasNoProd; hasNoSubs];

switch reactionType
    case {'both','all'} % For legacy reasons, 'both' is also allowed
        exchangeRxnsIndexes = allExch;
    case 'in'
        exchangeRxnsIndexes = hasNoSubs;
    case 'out'
        exchangeRxnsIndexes = hasNoProd;
    case 'blocked'
        exchangeRxnsIndexes = allExch(model.lb(allExch) == 0 & model.ub(allExch) == 0);
    case 'reverse'
        exchangeRxnsIndexes = allExch(model.lb(allExch) < 0 & model.ub(allExch) > 0);
    case 'uptake'
        exchangeRxnsIndexes = allExch([(model.lb(hasNoSubs) >= 0 & model.ub(hasNoSubs) > 0); ...
                              (model.lb(hasNoProd) < 0 & model.ub(hasNoProd) <= 0)]);
    case 'excrete'
        exchangeRxnsIndexes = allExch([(model.lb(hasNoSubs) < 0 & model.ub(hasNoSubs) <= 0); ...
                              (model.lb(hasNoProd) >= 0 & model.ub(hasNoProd) > 0)]);
    otherwise
        error('Invalid reactionType specified')
end
exchangeRxns = model.rxns(exchangeRxnsIndexes);
end
