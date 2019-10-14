% tutorial1
%   This script contains the code necessary for Exercise 1: importing a
%   model from Excel into a RAVEN model structure, as well as for running a
%   simulation with the parameters set in the Excel file.
%
%	Simonas Marcisauskas, 2019-10-14
%

%Import the Excel model into a RAVEN model structure
smallModel=importExcelModel('empty.xlsx');

%This solves the linear programming problem.
%NOTE: if sol.f is equal to zero then the problem is not solvable. Ensure
%that uptake and excretion of all necessary stuff are added to the model
%and run the solveLP again.
sol=solveLP(smallModel);

%Print the resulting exchange fluxes
printFluxes(smallModel,sol.x,true);

%Print all fluxes
printFluxes(smallModel,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
