function [growthRates, shadowPrices1, shadowPrices2] = runPhenotypePhasePlane(model, controlRxn1, controlRxn2, nPts, range1, range2)
% runPhenotypePhasePlane
%   Runs phenotype phase plane analysis and plots the results. The first
%   plot is a 3D surface plot showing the phenotype phase plane, the other
%   two plots show the shadow prices of the metabolites from the two
%   control reactions, which define the phases. Modified from the COBRA
%   phenotypePhasePlane function.
%
% Input:
%   model           a model structure
%   controlRxn1     reaction identifier of the first reaction to be plotted
%   controlRxn2     reaction identifier of the second reaction to be plotted
%   nPts            the number of points to plot in each dimension (opt,
%                   default 50)
%   range1          the range [from 0 to range1] of reaction 1 to plot
%                   (opt, default 20)
%   range2          the range [from 0 to range2] of reaction 2 to plot
%                   (opt, default 20)
%
% Output:
%   growthRates1    a matrix of maximum growth rates
%   shadowPrices1   a matrix with shadow prices for reaction 1
%   shadowPrices2   a matrix with shadow prices for reaction 2
%
% Modified from COBRA Toolbox phenotypePhasePlane.m
%
% Usage: [growthRates, shadowPrices1, shadowPrices2] = runPhenotypePhasePlane(model, controlRxn1, controlRxn2, nPts, range1, range2)
close all force % Close all existing figure windows (if open)
if nargin < 4
    nPts = 50;
end
if nargin < 5
    range1 = 20;
end
if nargin < 6
    range2 = 20;
end

rxnID1 = getIndexes(model,controlRxn1,'rxns',true);
metID1 = find(model.S(:,rxnID1));
rxnID2 = getIndexes(model,controlRxn2,'rxns',true);
metID2 = find(model.S(:,rxnID2));

% Create empty vectors for the results
ind1 = linspace(0,range1,nPts);
ind2 = linspace(0,range2,nPts);
growthRates = zeros(nPts);
shadowPrices1 = zeros(nPts);
shadowPrices2 = zeros(nPts);
[~,hsSol] = solveLP(model);
% Calulate points
PB = ProgressBar2(nPts,'Running PhPP analysis','cli');
for i = 1:nPts %ind1
    for j = 1:nPts %ind2
        model1 = setParam(model,'eq',controlRxn1,-1*ind1(i));
        model1 = setParam(model1,'eq',controlRxn2,-1*ind2(j));
        [fbasol,hsSol] = solveLP(model1,0,[],hsSol);
        try
            growthRates(j,i) = fbasol.x(logical(model1.c));
            shadowPrices1(j,i) = fbasol.sPrice(metID1);
            shadowPrices2(j,i) = fbasol.sPrice(metID2);   
        end
    end
    count(PB);
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');

% Plot the points
figure(2);
pcolor(ind1,ind2,shadowPrices1);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/hr)');
title(['Shadow price ' strrep(model.mets{metID1},'_','-')]);
colorbar();
xticklabels(sprintfc('%d', -xticks))
yticklabels(sprintfc('%d', -yticks))
figure(3);
pcolor(ind1,ind2,shadowPrices2);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/hr)');
title(['Shadow price ' strrep(model.mets{metID2},'_','-')]);
colorbar();
xticklabels(sprintfc('%d', -xticks))
yticklabels(sprintfc('%d', -yticks))
figure(1);
surfl(ind1,ind2,growthRates);
xlabel(strrep(strcat(controlRxn1,' (mmol/g DW-hr)'),'_','\_')), ylabel(strrep(strcat(controlRxn2,' (mmol/g DW-hr)'),'_','\_')), zlabel('growth rate (1/hr)');
colormap(hsv(17));
xticklabels(sprintfc('%d', -xticks))
yticklabels(sprintfc('%d', -yticks))
end
