function [minFluxes, maxFluxes, exitFlags]=getAllowedBounds(model,rxns,runParallel)
% getAllowedBounds
%   Returns the minimal and maximal fluxes through each reaction.
%
% Input:
%   model           a model structure
%   rxns            either a cell array of reaction IDs, a logical vector
%                   with the same number of elements as reactions in the
%                   model, or a vector of reaction indexes (opt, default
%                   model.rxns)
%   runParallel     speed up calculations by parallel processing. This is
%                   not beneficial if allowed bounds are calculated for
%                   only a few reactions, as the overhead of parallel
%                   processing will take longer. It requires MATLAB
%                   Parallel Computing Toolbox. If this is not installed,
%                   the calculations will not be parallelized, regardless
%                   what is indicated as runParallel. (opt, default true)
%
% Output:
%   minFluxes       minimal allowed fluxes
%   maxFluxes       maximal allowed fluxes
%   exitFlags       exit flags for min/max for each of the reactions. True
%                   if it was possible to calculate a flux
%
%   NOTE: In cases where no solution can be calculated, NaN is returned.
%
% Usage: [minFluxes, maxFluxes, exitFlags] = getAllowedBounds(model, rxns, runParallel)

if nargin<2 || isempty(rxns)
    rxns = 1:numel(model.rxns);
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns = convertCharArray(rxns);
    rxns = getIndexes(model,rxns, 'rxns');
end
if nargin<3
    runParallel = true;
end

[ps, oldPoolAutoCreateSetting] = parallelPoolRAVEN(runParallel);

minFluxes = zeros(numel(rxns),1);
maxFluxes = zeros(numel(rxns),1);
exitFlags = zeros(numel(rxns),2);
c = zeros(numel(model.rxns),1);

PB = ProgressBar2(numel(rxns),'Running getAllowedBounds','cli');
parfor i = 1:numel(rxns)
    count(PB)
    tmpModel = model;
    tmpModel.c = c;

    % Get minimal flux
    tmpModel.c(rxns(i)) = -1;
    solMin = solveLP(tmpModel);
    if ~isempty(solMin.f)
        minFluxes(i) = solMin.x(rxns(i));
    else
        minFluxes(i) = NaN;
    end

    % Get maximal flux
    tmpModel.c(rxns(i)) = 1;
    solMax=solveLP(tmpModel);
    exitFlags(i,:) = [solMin.stat solMax.stat];
    if ~isempty(solMax.f)
        maxFluxes(i) = solMax.x(rxns(i));
    else
        maxFluxes(i) = NaN;
    end
end
% Reset original Parallel setting
ps.Pool.AutoCreate = oldPoolAutoCreateSetting;
end
