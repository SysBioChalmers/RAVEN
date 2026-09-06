% tutorial3_solutions
%   This script contains the solutions for Tutorial 3, see Tutorial 3 on
%   the RAVEN wiki for more details:
%   https://github.com/SysBioChalmers/RAVEN/wiki/Tutorials
%   All the parameters are set in this script, rather than modifying the
%   Excel model file.

%Import the model
model=readYAMLmodel('smallYeast.yml');

%Step 1
%Set the upper bound of glucose uptake to 1 and O2 uptake to zero
model=setParam(model,'ub',{'glcIN' 'o2IN'},[1 0]);

%Set the objective to be ATP hydrolysis
model=setParam(model,'obj',{'ATPX'},1);

%Solve the model
sol=solveLP(model);

%Print the resulting fluxes. The ATP production rate should be 2.0. It was
%4.0 in Tutorial 2, but there sucrose was used instead of glucose.
printFluxes(model,sol.x,false);

%Step 2
%Check the yield of different products and print the results
%Change to fully aerobic
model=setParam(model,'ub',{'glcIN' 'o2IN'},[1 1000]);
model=setParam(model,'obj',{'ethOUT'},1);
sol=solveLP(model);
fprintf(['Yield of ethanol is ' num2str(sol.f) ' mol/mol\n']);
model=setParam(model,'obj',{'acOUT'},1);
sol=solveLP(model);
fprintf(['Yield of acetate is ' num2str(sol.f) ' mol/mol\n']);
model=setParam(model,'obj',{'glyOUT'},1);
sol=solveLP(model);
fprintf(['Yield of glycerol is ' num2str(sol.f) ' mol/mol\n']);
model=setParam(model,'obj',{'biomassOUT'},1);
sol=solveLP(model);
fprintf(['Yield of biomass is ' num2str(sol.f) '/h\n']);

%Step 3
%Solve for both aerobic and anaerobic growth
solA=solveLP(model);
model=setParam(model,'ub',{'o2IN'},0.5);
solB=solveLP(model);

%Step 4
%Change to anaerobic growth and maximize for biomass
model=setParam(model,'eq',{'o2IN'},0);
model=setParam(model,'obj',{'biomassOUT'},1);
sol=solveLP(model);
printFluxes(model,sol.x,true);
%One can see that the model predicts a glycerol production of 0.23
%mmol/gDW/h

%Run a single gene deletion
[genes, fluxes, originalGenes, details]=findGeneDeletions(model,'sgd');

%Get the indexes of these reactions
I=getIndexes(model,{'biomassOUT'},'rxns');
J=getIndexes(model,{'glyOUT'},'rxns');

okSolutions=find(fluxes(I,:)>10^-2); %Only look at solutions which are still growing
[maxGlycerol, J]=max(fluxes(J,okSolutions));
fprintf(['Glycerol production is ' num2str(maxGlycerol) ' after deletion of ' originalGenes{genes(okSolutions(J),:)} '\n']);

%The best gene deletion corresponds to turning off the ZWF1 reaction
%(YNL241C)
model2=setParam(model,'eq',{'ZWF'},0);
sol2=solveLP(model2);
compareFluxes(model,sol.x,sol2.x,'cutoff',10^-2, ...
    'metaboliteList',{'NADPH' 'NADH' 'NAD' 'NADP'});

%Step 5
%Read microarray results and calculate reporter metabolites (metabolites
%around which there are significant transcriptional changes)
[orfs, pvalues]=textread('expression.txt','%s%f');
repMets=reporterMetabolites(model,orfs,pvalues);
[I, J]=sort(repMets.metPValues);

fprintf('TOP 10 REPORTER METABOLITES:\n');
for i=1:min(numel(J),10)
    fprintf([repMets.mets{J(i)} '\t' num2str(I(i)) '\n']);
end
