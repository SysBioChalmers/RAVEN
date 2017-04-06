%This contains the code necessary for running Exercises 2. It is assumed
%that you have completed Exercises 1 before.

%Loads the model
model=importExcelModel('smallYeast.xlsx',true);

%Sets the upper bound of glucose uptake to 1 and O2 uptake to unlimited.
model=setParam(model,'ub',{'glcIN' 'o2IN'},[1 1000]);
%Set the objective to be ethanol production
model=setParam(model,'obj',{'ethOUT'},1);

%Solve the model
sol=solveLP(model);

%Print the resulting exchange fluxes
printFluxes(model,sol.x,true);

%This is for comparing two flux distributions
%Load the map
load 'pathway.mat' pathway;
drawMap('Aerobic vs Anaerobic',pathway,model,solA.x,solB.x,[],'mapFBA.pdf',10^-5);

%This is for running a single gene deletion
[genes fluxes originalGenes details]=findGeneDeletions(model,'sgd','fba');

%Get the indexes of these reactions
I=getIndexes(model,{'biomassOUT'},'rxns');
J=getIndexes(model,{'glyOUT'},'rxns');

okSolutions=find(fluxes(I,:)>10^-2); %Only look at solutions which are still growing
[maxGlycerol J]=max(fluxes(J,okSolutions));
maxGlycerol
originalGenes(genes(okSolutions(J),:))

%Draw map for the ZWF1 deletion strain
model2=setParam(model,'eq',{'ZWF'},0);
sol2=solveLP(model2);
drawMap('ZWF1 deletion vs WT',pathway,model,sol.x,sol2.x,[],'mapZWF.pdf',10^-5);
followChanged(model,sol2.x,sol.x, 10, 10^-2, 0,{'NADPH' 'NADH' 'NAD' 'NADP'});

%Load the model
SBMLFromExcel('smallYeast.xlsx','smallYeast.xml')
model=importModel('smallYeast.xml',true);
sol=solveLP(model);

%Define another model where all exchange reactions are open.
model2=model;
I=getIndexes(model,getExchangeRxns(model),'rxns');
model2.lb(I)=0;
model2.ub(I)=1000;

%Then delete ZWF
model2=setParam(model2,'eq',{'ZWF'},0);

%Then run MOMA
[fluxA,fluxB, flag]=qMOMA(model,model2);
drawMap('Aerobic vs Anaerobic MOMA',pathway,model,fluxA,fluxB,[],'mapMOMA.pdf',10^-5);

%Read microarray results and calculate reporter metabolites (metabolites
%around which there are significant transcriptional changes)
[orfs,pvalues]=textread('expression.txt','%s%f');
repMets=reporterMetabolites(model,orfs,pvalues);
[I J]=sort(repMets.metPValues);
 
fprintf('TOP 10 REPORTER METABOLITES:\n');
for i=1:min(numel(J),10)
   fprintf([repMets.mets{J(i)} '\t' num2str(I(i)) '\n']) 
end
 
%Get all reactions involving those metabolites and display them on a map
mets=ismember(model.mets,repMets.mets(J(1:10)));
[crap I]=find(model.S(mets,:));
pathway=trimPathway(pathway, model.rxns(I), true);
drawMap('Reactions involving the top 10 Reporter Metabolites',pathway,model,ones(numel(model.rxns),1),zeros(numel(model.rxns),1),[],'mapRM.pdf',10^-5);
