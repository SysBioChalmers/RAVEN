function [I, rxnNames]=getRxnsInComp(model,comp,includePartial)
% getRxnsInComp
%   Gets the reactions in a specified compartment
%
%   model           a model structure
%   comp            string with the compartment id
%   includePartial  true if reactions with metabolites in several
%                   compartments (normally transport reactions) should
%                   be included (optional, default false)
%
%   I               boolean vector of the reactions
%   rxnNames        the names of the reactions
%
% Usage: [I, rxnNames]=getRxnsInComp(model,comp,includePartial)

comp=char(comp);
if nargin<3
    includePartial=false;
end

J=find(ismember(upper(model.comps),upper(comp)));

if numel(J)~=1
    EM=['No unique match to compartment "' comp{1} '"'];
    dispEM(EM);
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
