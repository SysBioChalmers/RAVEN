function out=testLP(param)
modLoad('~all RAVENgit');
modelPath='/Users/hedani/Documents/Models/HMR/';

model=importModel([modelPath 'SCO_revised.xml'],false);
model=setParam(model,'obj',{'Biomass_SCO'},1);

setRavenSolver('gurobi');
gurobiRes=solveLP(model,1);
setRavenSolver('mosek');
mosekRes=solveLP(model,1);

out=[abs(mosekRes.x-gurobiRes.x)<0.01,mosekRes.x,gurobiRes.x];
array2table(out,'VariableNames',{'model' 'mosekRes' 'gurobiRes'})
end
