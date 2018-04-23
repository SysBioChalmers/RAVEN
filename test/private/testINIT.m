function out=testINIT(param)
modLoad('~all COBRA RAVENgit');
modelPath='/Users/hedani/Documents/Models/HMR/';

defParam=struct('solver','gurobi');
param=structUpdate(defParam,param);
setRavenSolver(param.solver);
changeCobraSolverParams('MILP','intTol',param.intTol);
changeCobraSolverParams('MILP','relMipGapTol',param.relMipGapTol);
changeCobraSolverParams('MILP','printLevel',param.printLevel);

%cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-4,'relMipGapTol',.01);
%cparams=structUpdate(cparams,params);

out=[];
hpaData=parseHPA([modelPath 'normal_tissue.csv']);
model=importModel([modelPath 'HMRdatabase2_00.xml'],false);

out=getINITModel(model,'colon',[],hpaData);
%colonModel=getINITModel(model,'colon',[],hpaData,[],[],'Dataset7.xlsx')
end

function s_merged=structUpdate(s_old,s_new)
s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
names = [fieldnames(s_merged); fieldnames(s_new)];
s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end
