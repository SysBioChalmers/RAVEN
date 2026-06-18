function [I, rxnNames]=getRxnsInComp(model,comp,varargin)
% getRxnsInComp  Get the reactions in a specified compartment.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% comp : char
%     string with the compartment id.
%
% Name-Value Arguments
% --------------------
% includePartial : logical
%     true if reactions with metabolites in several compartments (normally
%     transport reactions) should be included (default false).
%
% Returns
% -------
% I : double
%     boolean vector of the reactions.
% rxnNames : cell
%     the names of the reactions.
%
% Examples
% --------
%     [I, rxnNames] = getRxnsInComp(model, comp, includePartial);

comp=char(comp);
p=parseRAVENargs(varargin, {'includePartial',false});
includePartial=p.includePartial;

J=find(ismember(upper(model.comps),upper(comp)));

if numel(J)~=1
    EM=['No unique match to compartment "' comp{1} '"'];
    error('RAVEN:badInput', '%s', EM);
end

K=model.metComps==J; %Get all metabolites in the compartment

S=model.S~=0;

%Find the reactions which involve any of the mets
[~, I]=find(S(K,:));
I=unique(I);

%Then remove the ones which also include metabolites in other comps
if includePartial==false
    I=I(sum(S(:,I))==sum(S(K,I)));
end
rxnNames=model.rxnNames(I);
end
