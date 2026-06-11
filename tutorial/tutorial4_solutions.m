% tutorial4_solutions
%   This script contains the solutions for Tutorial 4, see Tutorial 4 in
%   "RAVEN tutorials.docx" for more details.
%
%   NOTE: Many of these changes are easier to do in the Excel sheet. They
%   are done here in code just to avoid having several model files.

%Import the model
model=readYAMLmodel('smallYeastBad.yml');

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

%Lots of ethanol produced. Also plot the equations to make the error easier
%to find
printFluxes(model,sol.x,false,[],[],'%rxnID (%rxnName):%flux\n\t%eqn\n');

%See that ADH1 should only produce one unit of ethanol. Change the reaction
%equation
model=changeRxns(model,'ADH1','acetaldehyde[c] + NADH[c] => ethanol[c] + NAD(+)[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,true); %Nothing produced, good

%Add excretion of all metabolites
model.b=[model.b inf(numel(model.b),1)];
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');

%By looking at the reactions which were unbalanced and that were in the
%flux list one can see that FBP should be changed to result in only one unit
%of F6P
model=changeRxns(model,'FBP','beta-D-fructofuranose 1,6-bisphosphate[c] => beta-D-fructofuranose 6-phosphate[c] + phosphate[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n');

%The same thing again and one should change PFK to only give one unit of
%F16P
model=changeRxns(model,'PFK','ATP[c] + beta-D-fructofuranose 6-phosphate[c] => ADP[c] + beta-D-fructofuranose 1,6-bisphosphate[c]',3);
sol=solveLP(model,1);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n'); %Now it works

%Set all uptakes and production to 0
model=setParam(model,'eq',getExchangeRxns(model),0);

%Since it is checked which metabolites could be consumed without
%production, one can no longer have free production of all metabolites
model.b=model.b(:,1);
I=canConsume(model);
disp(model.mets(I)); %These 12 metabolites can be consumed without any production

%Allow all uptake
model.b=[ones(numel(model.b),1)*-1000 model.b];

%Pick CO2 and force uptake of it
model=setParam(model,'eq',{'co2OUT'},-1); %Negative output means input
sol=solveLP(model);
printFluxes(model,sol.x,false,10^-5,[],'%rxnID (%rxnName):\n\t%eqn\n\t%flux\n'); %Now it works

%See that PDC converts pyruvate (3 carbons) to acetaldehyde (2 carbons)
%without any other products. If one googles, one may realize that CO2 is
%missing. This would be simpler to change in the Excel file (or using
%changeRxns), but one can change it here as an exercise. One therefore
%needs to find the index of the reactions and the index of cytosolic CO2 in
%order to change the reaction
Irxn=ismember(model.rxns,'PDC');
Imet=ismember(model.mets,'CO2_c');
model.S(Imet,Irxn)=1; %The coefficient is 1.0

%Display the new equation just to be sure
constructEquations(model,Irxn)

%The solution is now not feasible, meaning that it is no longer possible to
%force uptake of CO2 without any output
sol=solveLP(model);
