%This contains the code necessary for running Exercises 3. It is assumed
%that you are somewhat familiar with linear programming.

%NOTE: Many of these changes are easier to do in the Excel sheet. They are
%done here in code just to avoid having several model files.

%Load the model
model=importExcelModel('smallYeastBad.xlsx');

%Close all uptake and maximize for production
model=setParam(model,'eq',{'glcIN', 'o2IN'},[0 0]);
model=setParam(model,'obj',{'acOUT' 'biomassOUT' 'co2OUT' 'ethOUT' 'glyOUT'},[1 1 1 1 1]);
sol=solveLP(model);
printFluxes(model,sol.x,true); %Nothing produced, good

%Add some fake reactions
rxns.rxns={'FREE_ATP';'FREE_NADH';'FREE_NADPH'};
rxns.equations={'ATP <=> ADP + phosphate';'NAD(+) <=> NADH';'NADP(+) <=> NADPH'};
model=addRxns(model,rxns,2,'c');
sol=solveLP(model,1);
%Lots of ethanol produced.We also plot the equations to make the error easier to find
printFluxes(model,sol.x,false,[],[],'%rxnID (%rxnName):%flux\n\t%eqn\n');

%You will see that ADH1 should only produce one unit of ethanol. Change the
%reaction equation
model=changeRxns(model,'ADH1','acetaldehyde[c] + NADH[c] => ethanol[c] + NAD(+)[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,true); %Nothing produced, good

%Add excretion of all metabolites
model.b=[model.b inf(numel(model.b),1)];
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');

%By looking at the reactions which were unbalanced and that were in the
%flux list we can see that FBP should be changed to result in only one unit
%of F6P
model=changeRxns(model,'FBP','beta-D-fructofuranose 1,6-bisphosphate[c] => beta-D-fructofuranose 6-phosphate[c] + phosphate[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');
%The same thing again and we should change PFK to only give one unit of
%F16P
model=changeRxns(model,'PFK','ATP[c] + beta-D-fructofuranose 6-phosphate[c] => ADP[c] + beta-D-fructofuranose 1,6-bisphosphate[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n'); %Now it works

%Set all uptakes and production to 0
model=setParam(model,'eq',getExchangeRxns(model),0);

%Since we're checking for which metabolites could be consumed without
%production, we can no longer have free production of all metabolites.
model.b=model.b(:,1);
I=canConsume(model);
model.mets(I) %These 12 metabolites can be consumed without any production

%Allow all uptake
model.b=[ones(numel(model.b),1)*-1000 model.b];

%We pick CO2 and force uptake of it
model=setParam(model,'eq',{'co2OUT'},-1); %Negative output means input
sol=solveLP(model);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n'); %Now it works
%We see that PDC converts pyruvate (3 carbons) to acetaldehyde(2 carbons)
%without any other products. If we google we realize that CO2 is missing.
%This would be simpler to change in the Excel file (or using changeRxns),
%but lets change it here as an exercise. We need to find the index of the
%reactions and the index of cytosolic CO2 in order to change the reaction
Irxn=ismember(model.rxns,'PDC');
Imet=ismember(model.mets,'CO2_c');
model.S(Imet,Irxn)=1; %The coefficient is 1.0

%Display the new equation just to be sure
constructEquations(model,Irxn)

%The solution is now not feasible, meaning that it's no longer possible to
%force uptake of CO2 without any output
sol=solveLP(model)

%***Second part of tutorial
model=importExcelModel('smallYeastBad2.xlsx',true,false,true); %This has to be loaded with the setting to ignore error or it would find the error for you :)
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);
deletedReactions
deletedMetabolites

%It turned out that G15L_c was spelled G15Lc in one reaction. The best
%solution would be just to change it in the reaction list and remove the
%duplicate metabolite, but we can do it here as an exercise. We need the
%indexes of the two metabolites
Igood=ismember(model.mets,'G15L_c');
Ibad=ismember(model.mets,'G15Lc');

%Then get all reactions and the coefficients in which the wrong one participates
%move them to be for the right one instead
model.S(Igood,:)=model.S(Igood,:)+model.S(Ibad,:);

%Then delete the bad one
model=removeMets(model,'G15Lc');
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);
deletedReactions
deletedMetabolites

%The only difference was that we had 20 deleted metabolites instead of 21.
%Nothing too spectacular. Check production can tell us what we need to
%connect.
[notProducedMets, crap, neededForProductionMat,minToConnect]=checkProduction(model,true,model.comps,false);

%In order to have production of all 54 metabolites we need to enable
%production of these 12. This small model does not include net synthesis of
%co-factors, so let's concentrate on the other ones. Glycerone phosphate
%allows for connection 18 others, so it seems like a good target.
minToConnect

%If you google around a little bit, and know your metabolism, you would
%find that DHAP (dihydroxyacetone) and GLYP (glycerone phosphate) are
%actually synonymes. Lets only use DHAP
Igood=ismember(model.mets,'DHAP_c');
Ibad=ismember(model.mets,'GLYP_c');

%Then get all reactions and the coefficients in which the wrong one participates
%move them to be for the right one instead
model.S(Igood,:)=model.S(Igood,:)+model.S(Ibad,:);

%Then delete the bad one
model=removeMets(model,'GLYP_c');
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);
deletedReactions
deletedMetabolites
[notProducedMets, crap, neededForProductionMat,minToConnect]=checkProduction(model,true,model.comps,false);
minToConnect

%Still quite a lot of gaps and no immediate way to fix it. We could try
%including reactions from a reference network and see if that helps. Let's
%use the small yeast model from tutorial 2
refModel=importExcelModel('smallYeast.xlsx');
[newConnected cannotConnect addedRxns newModel]=fillGaps(model,{refModel},false);
addedRxns
newConnected

%By including the ALD6 reaction from the reference model we were able to
%connect 21 reactions!
[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(newModel,false,false,false,true);
deletedMetabolites
deletedReactions

%TADA! All the model seems to be connected

%All this stuff could be done in a more automated manner as well
model=importExcelModel('smallYeastBad2.xlsx',true,false,true);
gapReport(model,{refModel});
