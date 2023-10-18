% tutorial2
%   This exercise deals with the a small glycolysis model in RAVEN
%   compatible Excel format and shows the most basic aspects of the
%   stoichiometric modelling. It is shown how to build a simple model from
%   scratch, set parameters and perform simple simulations.
%   See Tutorial 2 in "RAVEN tutorials.docx" for more details.

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
