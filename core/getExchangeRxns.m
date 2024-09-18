function [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,reactionType)
% getExchangeRxns
%   Retrieves the exchange reactions from a model
%
%   model               a model structure
%   reactionType        retrieve all reactions ('both'), only production
%                       ('out'), or only consumption ('in') (optional, default
%                       'both')
%
%   exchangeRxns        cell array with the IDs of the exchange reactions
%   exchangeRxnsIndexes vector with the indexes of the exchange reactions
%
%   Exchange reactions are defined as reactions which involve only products
%   or only reactants. If the unconstrained field is present, then that is
%   used instead.
%
% Usage: [exchangeRxns,exchangeRxnsIndexes]=getExchangeRxns(model,reactionType)

if nargin<2
    reactionType='both';
else
    reactionType=char(reactionType);
end

hasNoProducts=sparse(numel(model.rxns),1);
hasNoReactants=sparse(numel(model.rxns),1);

if isfield(model,'unconstrained')
    if strcmpi(reactionType,'both') || strcmpi(reactionType,'out')
        [~, I]=find(model.S(model.unconstrained~=0,:)>0);
        hasNoProducts(I)=true;
    end
    if strcmpi(reactionType,'both') || strcmpi(reactionType,'in')
        [~, I]=find(model.S(model.unconstrained~=0,:)<0);
        hasNoReactants(I)=true;
    end
else
    if strcmpi(reactionType,'both') || strcmpi(reactionType,'out')
        hasNoProducts=sum((model.S>0))==0;
    end
    if strcmpi(reactionType,'both') || strcmpi(reactionType,'in')
        hasNoReactants=sum((model.S<0))==0;
    end
end
exchangeRxnsIndexes=find(hasNoProducts(:) | hasNoReactants(:));
exchangeRxns=model.rxns(exchangeRxnsIndexes);
end
