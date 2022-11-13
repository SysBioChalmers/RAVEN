function [biomassValues, targetValues, lineHandle] = runProductionEnvelope(model, deletions, lineColor, targetRxn, biomassRxn, geneDelFlag, nPts)
% runProductionEnvelope
%   Calculates the byproduct secretion envelope
%
% Input:
%   model            a model structure
%   targetRxn        identifier of target metabolite production reaction
%   biomassRxn       identifier of biomass reaction
%   nPts             number of points in the plot (opt, default 20)
%
% Output:
%   biomassValues    Biomass values for plotting
%   targetValues     Target upper and lower bounds for plotting
%
% Modified from COBRA Toolbox productionEnvelope.m
%
% Usage: [biomassValues, targetValues] = runProductionEnvelope(model, targetRxn, biomassRxn, nPts)

if nargin < 4
    nPts = 20;
end

% Run FBA to get upper bound for biomass
model   = setParam(model,'obj',biomassRxn,1);
[solMax,hsSol]  = solveLP(model);
solMax  = solMax.x(logical(model.c));
model.c = -model.c;
solMin  = solveLP(model,0,[],hsSol);
solMin  = solMin.x(logical(model.c));

% Create biomass range vector
biomassValues = linspace(solMin,solMax,nPts);

fprintf('Creating production envelope...   0%% complete');
% Max/min for target production
model = setParam(model,'obj',targetRxn,1);
for i = 1:length(biomassValues)
    progress=pad(num2str(floor(i/numel(biomassValues)*100)),3,'left');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);       
    model1 = setParam(model,'eq',biomassRxn,biomassValues(i));
    sol = solveLP(model1,0,[],hsSol);
    if (sol.stat > 0)
        targetUpperBound(i) = sol.x(logical(model.c));
    else
        targetUpperBound(i) = NaN;
    end
    model1.c = -model1.c;
    sol = solveLP(model1,0,[],hsSol);
    if (sol.stat > 0)
        targetLowerBound(i) = sol.x(logical(model1.c));
    else
        targetLowerBound(i) = NaN;
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');

% Plot results
lineHandle=plot([biomassValues fliplr(biomassValues)],[targetUpperBound fliplr(targetLowerBound)],lineColor,'LineWidth',2);
axis tight;
ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
xlabel('Growth rate (1/h)');

biomassValues = biomassValues';
targetValues = [targetLowerBound' targetUpperBound'];
