function out=testIO(param)
	modLoad('~all COBRA RAVENgit',[]);
	modelPath='/Users/hedani/Documents/Models/HMR/';

	defParam=struct('solver','gurobi');
	param=structUpdate(defParam,param);
	setRavenSolver(param.solver);
	
	out=[];
	
	disp('Testing importModel on RAVEN model')
	model=importModel([modelPath 'HMRdatabase2_00.xml'],false);
	sol=solveLP(model);

	disp('Testing importModel on COBRA model')
	model=importModel([modelPath 'HMRdatabase2_00_Cobra.xml'],false);
	sol=solveLP(model);
end

function s_merged=structUpdate(s_old,s_new)
	s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
	names = [fieldnames(s_merged); fieldnames(s_new)];
	s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end