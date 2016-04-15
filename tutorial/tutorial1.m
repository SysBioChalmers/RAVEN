%This contains the code necessary for importing a model from Excel and into
%a RAVEN model structure, as well as for running a simulation with the
%parameters set in the Excel file.

%This loads the Excel model and converts it into a RAVEN model structure
smallModel=importExcelModel('empty.xlsx')

%This solves the problem. If it looks like:
%sol=
%   f: []
%...
%then the problem is not solvable. Be sure that you have added uptake and
%excretion of all necessary stuff.
sol=solveLP(smallModel)

%Print the resulting exchange fluxes
printFluxes(smallModel,sol.x,true);

%Print all fluxes
printFluxes(smallModel,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');