function [minFluxes, maxFluxes, exitFlags]=getAllowedBounds(model,rxns)
% getAllowedBounds
%   Returns the minimal and maximal fluxes through each reaction.
%
%   model         a model structure
%   rxns          either a cell array of reaction IDs, a logical vector with the
%                 same number of elements as reactions in the model, or a vector
%                 of reaction indexes (opt, default model.rxns)
%
%   minFluxes     minimal allowed fluxes
%   maxFluxes     maximal allowed fluxes
%   exitFlags     exit flags for min/max for each of the reactions. True if
%                 it was possible to calculate a flux
%
%   NOTE: In cases where no solution can be calculated, NaN is returned.
%
%   Usage: [minFluxes, maxFluxes, exitFlags]=getAllowedBounds(model,rxns)
%
%   Rasmus Agren, 2013-04-21
%

if nargin<2
    rxns=1:numel(model.rxns);
else
    rxns=getIndexes(model,rxns, 'rxns');
end

minFluxes=zeros(numel(rxns),1);
maxFluxes=zeros(numel(rxns),1);
exitFlags=zeros(numel(rxns),2);

c=zeros(numel(model.rxns),1);
hsSolMin=[];
hsSolMax=[];
for i=1:numel(rxns)
    model.c=c;
    
    %Get minimal flux
    model.c(rxns(i))=-1;
    [solution, hsSolMin]=solveLP(model,0,[],hsSolMin);
    exitFlags(i,1)=solution.stat;
    if ~isempty(solution.f)
        minFluxes(i)=solution.x(rxns(i));
    else
        minFluxes(i)=NaN;
    end
    
    %Get maximal flux
    model.c(rxns(i))=1;
    [solution, hsSolMax]=solveLP(model,0,[],hsSolMax);
    exitFlags(i,2)=solution.stat;
    if ~isempty(solution.f)
        maxFluxes(i)=solution.x(rxns(i));
    else
        maxFluxes(i)=NaN;
    end
end
end
