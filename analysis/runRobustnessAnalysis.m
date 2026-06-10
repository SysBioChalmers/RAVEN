function [controlFlux, objFlux] = runRobustnessAnalysis(model, controlRxn, nPoints, objRxn, plotRedCost)
% runRobustnessAnalysis  Perform robustness analysis for a reaction and objective.
%
% Performs robustness analysis for a reaction of interest and an objective
% of interest.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% controlRxn : char
%     reaction of interest whose value is to be controlled.
% nPoints : double, optional
%     number of points to show on plot (default 20).
% objRxn : char, optional
%     reaction identifier of objective to be maximized (default uses the
%     objective defined in the model).
% plotRedCost : logical, optional
%     whether reduced cost should also be plotted (default false).
%
% Returns
% -------
% controlFlux : double
%     flux values of the reaction of interest, ranging from its minimum to
%     its maximum value.
% objFlux : double
%     optimal values of objective reaction at each control reaction flux
%     value.
%
% Examples
% --------
%     [controlFlux, objFlux] = runRobustnessAnalysis(model, controlRxn, ...
%         nPoints, objRxn);
%
% Notes
% -----
% Modified from COBRA Toolbox robustnessAnalysis.m.

if nargin < 3
    nPoints = 20;
end
if nargin < 4
    baseModel = model;
else
    baseModel = setParam(model,'obj',objRxn,1);
end
if nargin < 5
    plotRedCost = false;
end

if any(ismember(model.rxns,controlRxn))
    controlRxnIdx = getIndexes(model,controlRxn,'rxns');
    tmpModel = setParam(model,'obj',controlRxnIdx,1);
    solMax = solveLP(tmpModel);
    solMax = solMax.x(logical(tmpModel.c));
    tmpModel.c = -tmpModel.c;
    solMin = solveLP(tmpModel);
    solMin = solMin.x(logical(tmpModel.c));
else
    error('Control reaction does not exist!');
end

objFlux = zeros(nPoints,1);
redCost = zeros(nPoints,1);
controlFlux = linspace(solMin,solMax,nPoints)';


PB = ProgressBar2(length(controlFlux),'Running robustness analysis','cli');
for i=1:length(controlFlux)
    modelControlled = setParam(baseModel,'eq',controlRxnIdx,controlFlux(i));
    solControlled = solveLP(modelControlled);
    objFlux(i) = solControlled.x(logical(modelControlled.c));
    redCost(i) = solControlled.rCost(controlRxnIdx);
    count(PB);
end

if plotRedCost
    yyaxis right
    plot(controlFlux,redCost)
    ylabel([strrep(controlRxn,'_','-') ' reduced cost']);
    yyaxis left
end
plot(controlFlux,objFlux)
xlabel(strrep(controlRxn,'_','-'));
ylabel('Objective');
end
