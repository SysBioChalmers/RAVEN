function [concentrationMatrix, excRxnNames, timeVec, biomassVec] = runDynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
% runDynamicFBA
%   Performs dynamic FBA simulation using the static optimization approach
%
% Input:
%   model               a model structure
%   substrateRxns       cell array with exchange reaction identifiers for
%                       substrates that are initially in the media, whose
%                       concentration may change (e.g. not h2o or co2)
%   initConcentrations  numeric initial concentrations of substrates
%                       (matching substrateRxns)
%   initBiomass         numeric initial biomass (must be non-zero)
%   timeStep            numeric time step size
%   nSteps              numeric maximum number of time steps
%   plotRxns            cell array with exchange reaction identifiers for
%                       substrates whose concentration should be plotted
%   exclUptakeRxns      cell array with exchange reaction identifiers for
%                       substrates whose concentration does not change
%                       (e.g. co2, o2, h2o, h)
%
% Output:
%   concentrationMatrix numeric matrix with extracellular metabolite
%                       concentrations
%   excRxnNames         cell array with exchange reaction identifiers that
%                       match the metabolites included in the
%                       concentrationMatrix
%   timeVec             numeric vector of time points
%   biomassVec          numeric vector with biomass concentrations
%
% If no initial concentration is given for a substrate that has an open
% uptake in the model (i.e. `model.lb < 0`) the concentration is assumed to
% be high enough to not be limiting. If the uptake rate for a nutrient is
% calculated to exceed the maximum uptake rate for that nutrient specified
% in the model and the max uptake rate specified is > 0, the maximum uptake
% rate specified in the model is used instead of the calculated uptake
% rate.
%
% Modified from COBRA Toolbox dynamicFBA.m
%
% Usage: [concentrationMatrix, excRxnNames, timeVec, biomassVec] = runDynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)

% Find exchange rxns
excRxnNames = getExchangeRxns(model);
excRxnNames(ismember(excRxnNames,exclUptakeRxns))=[];
excInd = getIndexes(model,excRxnNames,'rxns');
% Figure out if substrate reactions are correct
missingInd = find(~ismember(substrateRxns,excRxnNames));
if (~isempty(missingInd))
    for i = 1:length(missingInd)
        fprintf('%s\n',substrateRxns{missingInd(i)});
    end
    error('Invalid substrate uptake reaction!');
end

% Initialize concentrations
[~, substrateMatchInd] = ismember(substrateRxns,excRxnNames);
concentrations = zeros(length(excRxnNames),1);
concentrations(substrateMatchInd) = initConcentrations;

% Deal with reactions for which there are no initial concentrations
originalBound = -model.lb(excInd);
noInitConcentration = (concentrations == 0 & originalBound > 0);
concentrations(noInitConcentration) = 1000;

biomass = initBiomass;

% Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

% Make sure bounds are not higher than what are specified in the model
aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
model.lb(excInd) = -uptakeBound;

concentrationMatrix = concentrations;
biomassVec = biomass;
timeVec(1) = 0;
[~,hsSol]=solveLP(model,1);
fprintf('Running dynamic FBA analysis...   0%% complete');
for stepNo = 1:nSteps
    % Run FBA
    [sol,hsSol] = solveLP(model,1,[],hsSol);
    mu = -sol.f;
    if (sol.stat ~= 1 || mu == 0)
        fprintf('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass);
        break;
    end
    uptakeFlux = sol.x(excInd);
    biomass = biomass*exp(mu*timeStep);
    biomassVec(end+1) = biomass;

    % Update concentrations
    concentrations = concentrations - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    concentrations(concentrations <= 0) = 0;
    concentrationMatrix(:,end+1) = concentrations;

    % Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
    % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    % Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    % Revert to original bounds if the rate was too high
    uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound(abs(uptakeBound) < 1e-9) = 0;

    model.lb(excInd) = -uptakeBound;
    timeVec(stepNo+1) = stepNo*timeStep;

    progress=pad(num2str(floor(stepNo/nSteps*100)),3,'left');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');

selNonZero = any(concentrationMatrix>0,2);
concentrationMatrix = concentrationMatrix(selNonZero,:);
excRxnNames = excRxnNames(selNonZero);
selPlot = ismember(excRxnNames,plotRxns);

% Plot concentrations as a function of time
clf
subplot(1,2,1);
plot(timeVec,biomassVec);
axis tight
title('Biomass');
xlabel('Time (h)');
ylabel('Concentration (g/L)');
subplot(1,2,2);
plot(timeVec,concentrationMatrix(selPlot,:));
axis tight
title('Substrates and/or products');
xlabel('Time (h)');
ylabel('Concentration (mmol/L)');
legend(strrep(excRxnNames(selPlot),'EX_',''));
end
