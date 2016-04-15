%This contains the code necessary for running Exercises 3. It is assumed
%that you are somewhat familiar with linear programming.

%Load the model
model=importExcelModel('smallYeastBad.xlsx');
model.b=[model.b inf(numel(model.b),1)];
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
I=canConsume(model);
model.mets(I)
gapReport(model,{refModel});

model=importExcelModel('smallYeastBad2.xlsx');
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);
deletedReactions
deletedMetabolites
[notProducedMets, crap, neededForProductionMat,minToConnect]=checkProduction(model,true,model.comps,false);
minToConnect

refModel=importExcelModel('smallYeast.xlsx');
[newConnected cannotConnect addedRxns newModel]=fillGaps(model,{refModel},false);
addedRxns
newConnected
cannotConnect