function [exchangeRxns, exchangeRxnsIndexes, exchangedMets]=getExchangeRxns(model,varargin)
% getExchangeRxns  Retrieve the exchange reactions from a model.
%
% Exchange reactions are identified by having either no substrates or no
% products.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% reactionType : char
%     which exchange reactions should be returned (default 'all'):
%
%     - 'all' : all reactions, irrespective of reaction bounds
%     - 'uptake' : reactions with bounds that imply that only uptake is
%       allowed. Reaction direction, upper and lower bounds are all
%       considered
%     - 'excrete' : reactions with bounds that imply that only excretion is
%       allowed. Reaction direction, upper and lower bounds are all
%       considered
%     - 'reverse' : reactions with non-zero upper and lower bounds that
%       imply that both uptake and excretion are allowed
%     - 'blocked' : reactions that have zero upper and lower bounds, not
%       allowing any flux
%     - 'in' : reactions where the boundary metabolite is the substrate of
%       the reaction; a positive flux value would imply uptake, but
%       reaction bounds are not considered
%     - 'out' : reactions where the boundary metabolite is the product of
%       the reaction; a negative flux value would imply uptake, but
%       reaction bounds are not considered
%
% Returns
% -------
% exchangeRxns : cell
%     cell array with the IDs of the exchange reactions.
% exchangeRxnsIndexes : double
%     vector with the indexes of the exchange reactions.
% exchangedMets : double
%     vector with the metabolite index for each exchange reaction (one
%     metabolite per exchange reaction). Only computed when requested.
%
% Notes
% -----
% The union of 'in' and 'out' equals 'all'. Also, the union of 'uptake',
% 'excrete', 'reverse' and 'blocked' equals 'all'.
%
% Examples
% --------
%     [exchangeRxns, exchangeRxnsIndexes] = getExchangeRxns(model, reactionType);

p=parseRAVENargs(varargin, {'reactionType','all'});
reactionType=char(p.reactionType);

% Find exchange reactions. Both branches must produce a column of reaction
% indexes: every use below (concatenation, bound lookups, model.rxns
% indexing) reads them as indexes, not as a mask over the reactions.
if isfield(model, 'unconstrained') && any(model.unconstrained~=0)
    % With boundary metabolites present, an exchange reaction is one that
    % involves a boundary metabolite. Producing the boundary metabolite is
    % what leaves the reaction without a product once the boundary
    % metabolite is removed, which is the definition the branch below uses.
    boundaryMets = model.unconstrained~=0;
    [~, I]=find(model.S(boundaryMets,:)>0);
    hasNoProd = unique(I(:));
    [~, I]=find(model.S(boundaryMets,:)<0);
    hasNoSubs = unique(I(:));
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
        
        exchangeRxnsIndexes = allExch([(model.lb(hasNoProd) < 0 & model.ub(hasNoProd) <= 0); ...
                              (model.lb(hasNoSubs) >= 0 & model.ub(hasNoSubs) > 0)]);
    case 'excrete'
        exchangeRxnsIndexes = allExch([(model.lb(hasNoProd) >= 0 & model.ub(hasNoProd) > 0); ...
                              (model.lb(hasNoSubs) < 0 & model.ub(hasNoSubs) <= 0)]);
    otherwise
        error('Invalid reactionType specified')
end
% unique rather than sort: a reaction that involves two boundary metabolites,
% and an empty reaction (no substrates and no products), each appear in both
% halves of allExch and would otherwise be reported twice.
exchangeRxnsIndexes = unique(exchangeRxnsIndexes);
exchangeRxns = model.rxns(exchangeRxnsIndexes);

if nargout > 2
    nExch=numel(exchangeRxnsIndexes);
    exchangedMets=zeros(nExch,1);
    for ei=1:nExch
        k=find(model.S(:,exchangeRxnsIndexes(ei)),1);
        if ~isempty(k)
            exchangedMets(ei)=k;
        end
    end
end
end
