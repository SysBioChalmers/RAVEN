% tutorial3
%   This exercise shows how to run FBA simulations, gene deletion analysis,
%   and how one can use GEMs as a scaffold for interpreting microarray
%   data. A simplified model of yeast metabolism is used in this approach
%   as an example.
%   See Tutorial 3 on the RAVEN wiki for more details:
%   https://github.com/SysBioChalmers/RAVEN/wiki/Tutorials
%
%   It is assumed that the user has already completed Tutorial 2

%Import the model
model=readYAMLmodel('smallYeast.yml');

%Set the upper bound of glucose uptake to 1 and O2 uptake to unlimited
model=setParam(model,'ub',{'glcIN' 'o2IN'},[1 1000]);

%Set the objective to be ethanol production
model=setParam(model,'obj',{'ethOUT'},1);

%Solve the model
sol=solveLP(model);

%Print the resulting exchange fluxes
printFluxes(model,sol.x,true);

%Run a single gene deletion
[genes, fluxes, originalGenes, details]=findGeneDeletions(model,'sgd');

%Get the indexes of these reactions
I=getIndexes(model,{'biomassOUT'},'rxns');
J=getIndexes(model,{'glyOUT'},'rxns');

okSolutions=find(fluxes(I,:)>10^-2); %Only look at solutions which are still growing
[maxGlycerol, J]=max(fluxes(J,okSolutions));
disp(maxGlycerol);
disp(originalGenes(genes(okSolutions(J),:)));

%Compare the ZWF1 deletion strain to wild type, looking only at the
%reactions that involve the redox cofactors
model2=setParam(model,'eq',{'ZWF'},0);
sol2=solveLP(model2);
compareFluxes(model,sol.x,sol2.x,'cutoff',10^-2, ...
    'metaboliteList',{'NADPH' 'NADH' 'NAD' 'NADP'});

%Reload the model, since earlier steps modified its bounds and objective
model=readYAMLmodel('smallYeast.yml');

%Read microarray results and calculate reporter metabolites (metabolites
%around which there are significant transcriptional changes)
[orfs, pvalues]=textread('expression.txt','%s%f');
repMets=reporterMetabolites(model,orfs,pvalues);
[I, J]=sort(repMets.metPValues);

fprintf('TOP 10 REPORTER METABOLITES:\n');
for i=1:min(numel(J),10)
    fprintf([repMets.mets{J(i)} '\t' num2str(I(i)) '\n']);
end
