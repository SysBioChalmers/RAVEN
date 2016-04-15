function [I metNames]=getMetsInComp(model,comp)
% getMetsInComp
%   Gets the metabolites in a specified compartment
%
%   model       a model structure
%   comp        string with the compartment id
%
%   I           boolean vector of the metabolites
%   metNames    the names of the metabolites
%
%   Usage: [I metNames]=getMetsInComp(model,comp)
%
%   Rasmus Agren, 2013-08-01
%

if ischar(comp)
    comp={comp};
end

J=find(ismember(upper(model.comps),upper(comp)));

if numel(J)~=1
   dispEM(['No unique match to compartment "' comp{1} '"']); 
end

I=model.metComps==J;
metNames=model.metNames(I);
end
