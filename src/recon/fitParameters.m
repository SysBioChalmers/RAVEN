function [parameters, fitnessScore, exitFlag, newModel]=fitParameters(model,xRxns,xValues,rxnsToFit,valuesToFit,parameterPositions,fitToRatio,initialGuess,plotFitting)
% fitParameters
%   Fits parameters such as maintenance ATP by quadratic programming
%
%   model                 a model structure
%   xRxns                 cell array with the IDs of the reactions that will be fixed for each data point
%   xValues               matrix with the corresponding values for each
%                         xRxns (columns are reactions)
%   rxnsToFit             cell array with the IDs of reactions that will be fitted to
%   valuesToFit           matrix with the corresponding values for each
%                         rxnsToFit (columns are reactions)
%   parameterPositions    stucture that determines where the parameters are in the
%                         stoichiometric matrix. Contains the fields:
%   	position          cell array of vectors where each element contains
%                         the positions in the S-matrix for that parameter
%   	isNegative        cell array of vectors where the elements are true
%                         if that position should be the negative of the
%                         fitted value (to differentiate between
%                         production/consumption)
%	fitToRatio            if the ratio of simulated to measured values should
%                         be fitted instead of the absolute value. Used to prevent
%                         large fluxes from having too large impact (opt,
%                         default true)
%   initialGuess          initial guess of the parameters (opt)
%   plotFitting           true if the resulting fitting should be plotted
%                         (opt, default false)
%
%   parameters            fitted parameters in the same order as in
%                         parameterPositions
%   fitnessScore          the correponding residual sum of squares
%   newModel              updated model structure with the fitted parameters
%
%   Usage: [parameters, fitnessScore, exitFlag, newModel]=fitParameters(model,...
%           xRxns,xValues,rxnsToFit,valuesToFit,parameterPositions,fitToRatio,...
%           initialGuess,plotFitting)

if nargin<7
    fitToRatio=true;
end
if nargin<8
    initialGuess=ones(numel(parameterPositions.position),1);
end
if isempty(initialGuess)
    initialGuess=ones(numel(parameterPositions.position),1);
end
if nargin<9
    plotFitting=false;
end

%Find the indexes of reactions that will be fitted
[I, rxnsToFitIndexes]=ismember(rxnsToFit,model.rxns);

if ~all(I)
    EM='Could not find all reactions in rxnsToFit';
    dispEM(EM);
end

%Find the indexes of reactions that will be used for constraints.
[I, xRxnsIndexes]=ismember(xRxns,model.rxns);

if ~all(I)
    EM='Could not find all reactions in xRxns';
    dispEM(EM);
end

[parameters, fitnessScore, exitFlag]=fminsearch(@(parameters) getRSS(parameters,model,xRxnsIndexes,xValues,rxnsToFitIndexes,valuesToFit,parameterPositions,fitToRatio),initialGuess);

parameters=abs(parameters);

if plotFitting==true
    %Set the resulting parameters
    [~, resultingFluxes, newModel]=getRSS(parameters,model,xRxnsIndexes,xValues,rxnsToFitIndexes,valuesToFit,parameterPositions,true);
    plot(xValues,valuesToFit,'o',xValues,resultingFluxes,'-*');
end
end

function [rss, resultingFluxes, newModel]=getRSS(parameters,model,xRxnsIndexes,xValues,rxnsToFitIndexes,valuesToFit,parameterPositions,fitToRatio)
parameters=abs(parameters);

%Set the parameters at the positions specified in parameterPositions
for i=1:numel(parameterPositions.position)
    %Set positive
    model.S(parameterPositions.position{i}(parameterPositions.isNegative{i}==false))=parameters(i);
    
    %Set negative
    model.S(parameterPositions.position{i}(parameterPositions.isNegative{i}==true))=parameters(i)*-1;
end

%Also return an updated model
newModel=model;

%Loop through each data point, set xRxns to xValues and calculate the sum
%of squares for the rxnsToFit
rss=0;
resultingFluxes=[];
for i=1:size(xValues,1)
    %Fix for more xRxns!
    model.lb(xRxnsIndexes)=xValues(i,:);
    model.ub(xRxnsIndexes)=xValues(i);
    
    sol=solveLP(model);
    
    %Calculate the rss
    if fitToRatio==false
        rs=sol.x(rxnsToFitIndexes)'-valuesToFit(i,:);
    else
        rs=sol.x(rxnsToFitIndexes)'./valuesToFit(i,:)-ones(1,size(valuesToFit,2));
    end
    rss=rss+rs*rs';
    resultingFluxes=[resultingFluxes sol.x(rxnsToFitIndexes)];
end
end
