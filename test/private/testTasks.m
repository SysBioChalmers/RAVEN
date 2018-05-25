function out=testTasks(param)
	modLoad('~all COBRA RAVENgit',[]);
	modelPath='/Users/hedani/Documents/Models/HMR/';

	defParam=struct('solver','gurobi');
	param=structUpdate(defParam,param);
	setRavenSolver(param.solver);
	
	out=[];
	hpaData=parseHPA([modelPath 'normal_tissue.csv']);
	%model=importModel([modelPath 'colonHedani.xml'],false);
	model=importModel([modelPath 'HMRdatabase2_00.xml'],false);

	%[taskReport1 essentialRxns1 taskStructure1]=checkTasks(model,[modelPath 'Dataset7.xlsx'],true,false,true);

	modelInit=getINITModel(model,'colon',[],hpaData,[],[],[modelPath 'Dataset7.xlsx']);

	[taskReport2, essentialRxns2, taskStructure2]=checkTasks(modelInit,[modelPath 'Dataset7.xlsx'],true,false,true);
	out=modelInit;
end

function s_merged=structUpdate(s_old,s_new)
	s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
	names = [fieldnames(s_merged); fieldnames(s_new)];
	s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end
