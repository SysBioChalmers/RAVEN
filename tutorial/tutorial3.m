% tutorial3
%   This exercise shows how to run FBA and minimization of metabolic
%   adjustment (MOMA) simulations and how one can use GEMs as a scaffold
%   for interpreting microarray data. A simplified model of yeast
%   metabolism is used in this approach as an example.
%   See Tutorial 3 in "RAVEN tutorials.docx" for more details.
%
%   It is assumed that the user has already completed Tutorial 2

%Import the Excel model
model=importExcelModel('smallYeast.xlsx',true);

%Set the upper bound of glucose uptake to 1 and O2 uptake to unlimited
model=setParam(model,'ub',{'glcIN' 'o2IN'},[1 1000]);

%Set the objective to be ethanol production
model=setParam(model,'obj',{'ethOUT'},1);

%Solve the model
sol=solveLP(model);

%Print the resulting exchange fluxes
printFluxes(model,sol.x,true);

%Compare two flux distributions by loading the map
load 'pathway.mat' pathway;
drawMap('Aerobic vs Anaerobic',pathway,model,solA.x,solB.x,[],'mapFBA.pdf',10^-5);

%Run a single gene deletion
[genes, fluxes, originalGenes, details]=findGeneDeletions(model,'sgd','fba');

%Get the indexes of these reactions
I=getIndexes(model,{'biomassOUT'},'rxns');
J=getIndexes(model,{'glyOUT'},'rxns');

okSolutions=find(fluxes(I,:)>10^-2); %Only look at solutions which are still growing
[maxGlycerol, J]=max(fluxes(J,okSolutions));
disp(maxGlycerol);
disp(originalGenes(genes(okSolutions(J),:)));

%Draw map for the ZWF1 deletion strain
model2=setParam(model,'eq',{'ZWF'},0);
sol2=solveLP(model2);
drawMap('ZWF1 deletion vs WT',pathway,model,sol.x,sol2.x,[],'mapZWF.pdf',10^-5);
followChanged(model,sol2.x,sol.x, 10, 10^-2, 0,{'NADPH' 'NADH' 'NAD' 'NADP'});

%Import the model
SBMLFromExcel('smallYeast.xlsx','smallYeast.xml')
model=importModel('smallYeast.xml',true);
sol=solveLP(model);

%Define another model where all exchange reactions are open
model2=model;
I=getIndexes(model,getExchangeRxns(model),'rxns');
model2.lb(I)=0;
model2.ub(I)=1000;

%Delete ZWF gene
model2=setParam(model2,'eq',{'ZWF'},0);

%Run MOMA
[fluxA, fluxB, flag]=qMOMA(model,model2);
drawMap('Aerobic vs Anaerobic MOMA',pathway,model,fluxA,fluxB,[],'mapMOMA.pdf',10^-5);

%Read microarray results and calculate reporter metabolites (metabolites
%around which there are significant transcriptional changes)
[orfs, pvalues]=textread('expression.txt','%s%f');
repMets=reporterMetabolites(model,orfs,pvalues);
[I, J]=sort(repMets.metPValues);

fprintf('TOP 10 REPORTER METABOLITES:\n');
for i=1:min(numel(J),10)
    fprintf([repMets.mets{J(i)} '\t' num2str(I(i)) '\n']);
end

%Get all reactions involving those metabolites and display them on a map
mets=ismember(model.mets,repMets.mets(J(1:10)));
[~, I]=find(model.S(mets,:));
pathway=trimPathway(pathway, model.rxns(I), true);
drawMap('Reactions involving the top 10 Reporter Metabolites',pathway,model,ones(numel(model.rxns),1),zeros(numel(model.rxns),1),[],'mapRM.pdf',10^-5);
