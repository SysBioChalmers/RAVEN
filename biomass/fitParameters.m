function [parameters, fitnessScore, exitFlag, newModel]=fitParameters(model,xRxns,xValues,rxnsToFit,valuesToFit,parameterPositions,varargin)
% fitParameters  Fit parameters such as maintenance ATP by quadratic programming.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% xRxns : cell
%     cell array with the IDs of the reactions that will be fixed for each
%     data point.
% xValues : double
%     matrix with the corresponding values for each xRxns (columns are
%     reactions).
% rxnsToFit : cell
%     cell array with the IDs of reactions that will be fitted to.
% valuesToFit : double
%     matrix with the corresponding values for each rxnsToFit (columns are
%     reactions).
% parameterPositions : struct
%     structure that determines where the parameters are in the
%     stoichiometric matrix, with fields:
%
%     - position : cell array of vectors where each element contains the
%       positions in the S-matrix for that parameter
%     - isNegative : cell array of vectors where the elements are true if
%       that position should be the negative of the fitted value (to
%       differentiate between production/consumption)
%
% Name-Value Arguments
% --------------------
% fitToRatio : logical
%     if the ratio of simulated to measured values should be fitted
%     instead of the absolute value. Used to prevent large fluxes from
%     having too large an impact (default true).
% initialGuess : double
%     initial guess of the parameters (default ones).
% plotFitting : logical
%     true if the resulting fitting should be plotted (default false).
%
% Returns
% -------
% parameters : double
%     fitted parameters in the same order as in parameterPositions.
% fitnessScore : double
%     the corresponding residual sum of squares.
% exitFlag : double
%     exit status returned by fminsearch.
% newModel : struct
%     updated model structure with the fitted parameters.
%
% Examples
% --------
%     [parameters, fitnessScore, exitFlag, newModel]=fitParameters(model,...
%         xRxns,xValues,rxnsToFit,valuesToFit,parameterPositions,fitToRatio,...
%         initialGuess,plotFitting);

p=parseRAVENargs(varargin, {'fitToRatio',true; 'initialGuess',[]; 'plotFitting',false});
fitToRatio=p.fitToRatio;
initialGuess=p.initialGuess;
if isempty(initialGuess)
    initialGuess=ones(numel(parameterPositions.position),1);
end
plotFitting=p.plotFitting;

xRxns=convertCharArray(xRxns);
rxnsToFit=convertCharArray(rxnsToFit);

%Find the indexes of reactions that will be fitted
[I, rxnsToFitIndexes]=ismember(rxnsToFit,model.rxns);

if ~all(I)
    EM='Could not find all reactions in rxnsToFit';
    error('RAVEN:badInput', '%s', EM);
end

%Find the indexes of reactions that will be used for constraints.
[I, xRxnsIndexes]=ismember(xRxns,model.rxns);

if ~all(I)
    EM='Could not find all reactions in xRxns';
    error('RAVEN:badInput', '%s', EM);
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
