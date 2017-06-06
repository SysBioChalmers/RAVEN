%This is the solution to tutorial 1. The file small.xlsx contains all the
%necessary modifications to solve for ATP production in glycolysis. 
%
% Rasmus Agren, 2013-08-06
% Simonas Marcisauskas, 2017-06-06 - revision
%

%This loads the Excel model and converts it into a RAVEN model structure
smallModel=importExcelModel('small.xlsx');

%This solves the problem. If it looks like:
%sol=
%   f: []
%...
%then the problem is not solvable. Be sure that you have added uptake and
%excretion of all necessary stuff.
sol=solveLP(smallModel);

%Print the resulting exchange fluxes
printFluxes(smallModel,sol.x,true);

%Print all fluxes
printFluxes(smallModel,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
