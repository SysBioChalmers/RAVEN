function [fluxA,fluxB, flag]=qMOMA(modelA,modelB,fluxMinWeight)
% qMOMA
%   Uses quadratic programming to minimize the sum((fluxAi - fluxBi)^2)
%
%   modelA        a model structure for the test case. This model must be a
%                 subset of modelB (no reactions that are not in modelB)
%   modelB        a model structure for the reference case
%   fluxMinWeight a double >=1 that determines whether minimization of the
%                 sum of fluxes should also be taken into account in the
%                 optimization. A value of 2.0 means that sum(fluxAi)^2 +
%                 sum(fluxBi)^2 has equal weight as sum((fluxAi - fluxBi)^2).
%                 Values of around 1.01 should be enough to get rid of loops
%                 (opt, default 1)
%
%   fluxA         the resulting flux distribution in the test model
%   fluxB         the resulting flux distribution in the reference model
%   flag          1 if the optimization terminated successfully
%
%   Usage: [fluxA,fluxB, flag]=qMOMA(modelA,modelB,fluxMinWeight)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<3
    fluxMinWeight=1;
end

%Match the reactions and metabolites in the small model (modelA) to the
%larger model
[rxnExists,mapRxns]=ismember(modelA.rxns,modelB.rxns);

%Check that the smaller model is a subset of the larger one
if any(~rxnExists)
    dispEM('All reactions in the test model must exist in the reference model');
end

%In order to make the calculations a little easier to formulate I reshape
%modelA.S so that it is the same dimension as modelB.S
S=modelA.S;
lb=modelA.lb;
ub=modelA.ub;
c=modelA.c;
modelA.S=sparse(numel(modelA.mets),numel(modelB.rxns));
modelA.lb=zeros(numel(modelB.rxns),1);
modelA.ub=modelA.lb;
modelA.c=modelA.lb;
modelA.S(:,mapRxns)=S;
modelA.lb(mapRxns)=lb;
modelA.ub(mapRxns)=ub;
modelA.c(mapRxns)=c;

%Create the H matrix for the quadratic problem. Note that there are no
%linear terms in the equation

%Equation to solve is:(xA1 - xB1)^2 + (xA2 - xB2)^2 ....
%Is equal to xA1^2 - xA1*xB1 - xB1*xA1 + xB1^2 + xA2^2 - xA2*xB2 - xB2*xA2 + xB2^2...

%For three fluxes this would give the following H matrix (first three rows
%for A and last three rows for B)

%   A1 A2 A3 B1 B2 B3
%A1 1  0  0  -1  0  0
%A2 0  1  0  0  -1  0
%A3 0  0  1  0  0  -1
%B1 -1  0  0  1  0  0
%B2 0  -1  0  0  1  0
%B3 0  0  -1  0  0  1

%Create stoichiometric matrix and bounds that contain both models
fullS=[modelA.S sparse(size(modelA.S,1),size(modelA.S,2));sparse(size(modelB.S,1),size(modelB.S,2)) modelB.S];
fullLB=[modelA.lb;modelB.lb];
fullUB=[modelA.ub;modelB.ub];
fullB=zeros(size(modelA.S,1)+size(modelB.S,1),1);

H=[eye(size(fullS,2)/2)*fluxMinWeight eye(size(fullS,2)/2)*-1;eye(size(fullS,2)/2)*-1 eye(size(fullS,2)/2)*fluxMinWeight];

x=quadprog(H,zeros(size(H,1),1),[],[],fullS,fullB,fullLB,fullUB);

if any(x)
    fluxA=x(1:numel(modelB.rxns));
    fluxA=fluxA(mapRxns);
    fluxB=x(numel(modelB.rxns)+1:end);
    flag=1; %Since it never converges good enough
else
    fluxA=zeros(numel(modelA.rxns),1); %Still the old number
    fluxB=zeros(numel(modelB.rxns),1);
    flag=-1;
end
end
