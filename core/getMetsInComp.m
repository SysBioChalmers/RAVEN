function [I, metNames]=getMetsInComp(model,comp)
% getMetsInComp
%   Gets the metabolites in a specified compartment
%
%   model       a model structure
%   comp        string with the compartment id
%
%   I           boolean vector of the metabolites
%   metNames    the names of the metabolites
%
%   Usage: [I, metNames]=getMetsInComp(model,comp)

comp=char(comp);

J=find(ismember(upper(model.comps),upper(comp)));

if numel(J)~=1
    EM=['No unique match to compartment "' comp{1} '"'];
    dispEM(EM);
end

I=model.metComps==J;
metNames=model.metNames(I);
end
