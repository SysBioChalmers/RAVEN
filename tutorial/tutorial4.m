% tutorial4
%   This script contains the list of functions necessary for running
%   Tutorial 4, see Tutorial 4 in "RAVEN tutorials.docx" for more details.
%   Several key stages may be missing, try to fill these gaps before
%   checking the solutions in tutorial4_solutions. It is assumed that the
%   user is somewhat familiar with linear programming.
%
%   It is assumed that the user has already completed Tutorials 2-3

%Import the model to gap-fill and the reference model to draw reactions from
model=readYAMLmodel('smallYeastBad.yml');
refModel=readYAMLmodel('smallYeast.yml');
model.b=[model.b inf(numel(model.b),1)];
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
I=canConsume(model);
disp(model.mets(I));
gapReport(model,{refModel});

[newConnected, cannotConnect, addedRxns, newModel]=fillGaps(model,{refModel},false);
disp(addedRxns);
disp(newConnected);
disp(cannotConnect);
