% tutorial3
%   This script contains the list of functions necessary for running
%   Exercise 3, see Exercise 3 in "RAVEN tutorials.docx" for more
%   details. Several key stages may be missing, try to fill these gaps
%   before checking the solutions in tutorial3_solutions. It is assumed
%   that the user is somewhat familiar with linear programming.
%
%   It is assumed that the user has already completed Exercises 1-2
%   (tutorial1-tutorial2)
%
%   Simonas Marcisauskas, 2019-10-21
%

%Import the Excel model
model=importExcelModel('smallYeastBad.xlsx');
model.b=[model.b inf(numel(model.b),1)];
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
I=canConsume(model);
disp(model.mets(I));
gapReport(model,{refModel});

model=importExcelModel('smallYeastBad2.xlsx');
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);
disp(deletedReactions);
disp(deletedMetabolites);
[notProducedMets, ~, neededForProductionMat, minToConnect]=checkProduction(model,true,model.comps,false);
disp(minToConnect);

refModel=importExcelModel('smallYeast.xlsx');
[newConnected, cannotConnect, addedRxns, newModel]=fillGaps(model,{refModel},false);
disp(addedRxns);
disp(newConnected);
disp(cannotConnect);
