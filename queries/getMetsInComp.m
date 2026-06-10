function [I, metNames]=getMetsInComp(model,comp)
% getMetsInComp  Get the metabolites in a specified compartment.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% comp : char
%     string with the compartment id.
%
% Returns
% -------
% I : logical
%     boolean vector of the metabolites.
% metNames : cell
%     the names of the metabolites.
%
% Examples
% --------
%     [I, metNames] = getMetsInComp(model, comp);

comp=char(comp);

J=find(ismember(upper(model.comps),upper(comp)));

if numel(J)~=1
    EM=['No unique match to compartment "' comp{1} '"'];
    dispEM(EM);
end

I=model.metComps==J;
metNames=model.metNames(I);
end
