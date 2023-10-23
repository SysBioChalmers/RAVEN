function [minFluxes, maxFluxes, exitFlags]=getAllowedBounds(model,rxns,runParallel)
% getAllowedBounds
%   Returns the minimal and maximal fluxes through each reaction.
%
% Input:
%   model       a model structure
%   rxns        either a cell array of reaction IDs, a logical vector with the
%               same number of elements as reactions in the model, or a vector
%               of reaction indexes (opt, default model.rxns)
%   runParallel make use of MATLAB parallel pool to speed up calculations. Not
%               beneficial if only a limited number of reactions are simulated.
%               (opt, default true)
%
% Output:
%   minFluxes   minimal allowed fluxes
%   maxFluxes   maximal allowed fluxes
%   exitFlags   exit flags for min/max for each of the reactions. True if it was
%               possible to calculate a flux
%
%   NOTE: In cases where no solution can be calculated, NaN is returned.
%
% Usage: [minFluxes, maxFluxes, exitFlags] = getAllowedBounds(model, rxns, runParallel)

if nargin<2
    rxns = 1:numel(model.rxns);
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns = convertCharArray(rxns);
    rxns = getIndexes(model,rxns, 'rxns');
end
if nargin<3
    runParallel = true;
end
if runParallel
    addonList = matlab.addons.installedAddons;
    if ~any(strcmpi(addonList.Name,'Parallel Computing Toolbox'))
        disp('Cannot find MATLAB Parallel Computing Toolbox, process is not parallelized.')
        runParallel = false;
    end
end

minFluxes = zeros(numel(rxns),1);
maxFluxes = zeros(numel(rxns),1);
exitFlags = zeros(numel(rxns),2);
c = zeros(numel(model.rxns),1);

p = 1;
progressbar('Calculating the minimal and maximal fluxes')
if runParallel
    D = parallel.pool.DataQueue;
    afterEach(D, @updateProgressParallel);
    parfor i = 1:numel(rxns)
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
        send(D, i);
    end
else
    for i = 1:numel(rxns)
        tmpModel = model;
        tmpModel.c = c;

        % Get minimal flux
        tmpModel.c(rxns(i)) = -1;
        solMin=solveLP(tmpModel);
        if ~isempty(solMin.f)
            minFluxes(i) = solMin.x(rxns(i));
        else
            minFluxes(i) = NaN;
        end

        % Get maximal flux
        tmpModel.c(rxns(i)) = 1;
        solMax = solveLP(tmpModel);
        exitFlags(i,:) = [solMin.stat solMax.stat];
        if ~isempty(solMax.f)
            maxFluxes(i) = solMax.x(rxns(i));
        else
            maxFluxes(i) = NaN;
        end
        progressbar(i/numel(rxns))
    end
end
progressbar(1) % Make sure it closes

    function updateProgressParallel(~)
        progressbar(p/numel(rxns))
        p = p + 1;
    end
end
