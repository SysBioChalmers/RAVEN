function [controlFlux, objFlux] = runRobustnessAnalysis(model, controlRxn, nPoints, objRxn)
% runRobustnessAnalysis
%   Performs robustness analysis for a reaction of interest and an objective
%   of interest. Modified from the COBRA robustnessAnalysis function.
%
% Input:
%   model           a model structure
%   controlRxn      reaction of interest whose value is to be controlled
%   nPoints         number of points to show on plot (opt, default 20)
%   objRxn          reaction identifier of objective to be maximized (opt,
%                   default it uses the objective defined in the model)
%
% Output:
%   controlFlux     flux values of the reaction of interest, ranging from
%                   its minimum to its maximum value
%   objFlux         optimal values of objective reaction at each control
%                   reaction flux value
%
% Usage: runRobustnessAnalysis(model, controlRxn, nPoints, objRxn)

if nargin < 3
    nPoints = 20;
end
if nargin > 4
    baseModel = setParam(model,'obj',objRxn,1);
else
    baseModel = model;
end

if any(ismember(model.rxns,controlRxn))
    tmpModel = setParam(model,'obj',controlRxn,1);
    solMax = solveLP(tmpModel);
    solMax = solMax.x(logical(tmpModel.c));
    tmpModel.c = -tmpModel.c;
    solMin = solveLP(tmpModel);
    solMin = solMin.x(logical(tmpModel.c));
else
    error('Control reaction does not exist!');
end

objFlux = zeros(nPoints,1);
controlFlux = linspace(solMin,solMax,nPoints)';

fprintf('Running robustness analysis...   0%% complete');
for i=1:length(controlFlux)
    progress=pad(num2str(floor(i/numel(controlFlux)*100)),3,'left');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);    
    modelControlled = setParam(baseModel,'eq',controlRxn,controlFlux(i));
    solControlled = solveLP(modelControlled);
    objFlux(i) = solControlled.x(logical(modelControlled.c));
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');

plot(controlFlux,objFlux)
xlabel(strrep(controlRxn,'_','-'));
ylabel('Objective');
axis tight;
end
