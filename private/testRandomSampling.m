function out=testRandomSampling(param)
	%modLoad('~all COBRA RAVENgit')
	modelPath='/Users/hedani/Documents/Models/HMR/';

	model=importModel([modelPath 'SCO_revised.xml'],false);
	model=setParam(model,'obj',{'Biomass_SCO'},1);

	setRavenSolver('gurobi');
	[res crap prob]=solveLP(model);
	ts=tic;
	for i=1:100 solveLP(model); end
	tTot=toc(ts);


	gparams=struct('Presolve',2,'TimeLimit',1e9,'OutputFlag',0,'MIPGap',.01);
	%gparams=structUpdate(gparams,param);
	gprob=mosekToGurobiProb(prob); 
	ts=tic;
	for i=1:100 mosekToGurobiProb(prob); end
	tConv=toc(ts);

	res=gurobi(gprob,gparams);
	ts=tic;
	for i=1:100 gurobi(gprob,gparams); end
	tSolver=toc(ts);

	ts=tic;
	for i=1:100 gurobiToMosekRes(res,length(prob.c),false); end
	tConv=tConv+toc(ts);

	tOverHead=tTot-tConv-tSolver;

	out=[tSolver tConv tOverHead];

	setRavenSolver('mosek');
	params.relGap=0.4;
	params.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_FREE_SIMPLEX';

	ts=tic;
	for i=1:100 mosekopt(['minimize echo(0)'],prob,getMILPParams(params)); end
	tMskSolve=toc(ts);

	ts=tic;
	for i=1:100 solveLP(model); end
	tOverHead2=toc(ts)-tMskSolve;

	out=[out;tMskSolve 0 tOverHead2];

	array2table(out,'VariableNames',{'tSolver','tStructConv','tOther'})
	%profile on
	%mosekToGurobiProb(prob);
	%gurobiToMosekRes(res,length(prob.c),false);
	%profile viewer
end