function out=testLP2(param)
	modLoad('~all COBRA RAVENgit')
	modelPath='/Users/hedani/Documents/Models/HMR/';

	model=importModel([modelPath 'SCO_revised.xml'],false);
	model=setParam(model,'obj',{'Biomass_SCO'},1);

	setRavenSolver('gurobi');
	gurobiRes=solveLP(model,1);
	setRavenSolver('mosek');
	mosekRes=solveLP(model,1);
	setRavenSolver('cobra');
	cobraRes=solveLP(model,1);

	out.f=[mosekRes.f,gurobiRes.f,cobraRes.f];
	out.similarity=[arrSimilarity(cobraRes.x,gurobiRes.x,0.01) arrSimilarity(mosekRes.x,gurobiRes.x,0.01) arrSimilarity(mosekRes.x,cobraRes.x,0.01)];
	out.x=[mosekRes.x,gurobiRes.x,cobraRes.x];
	array2table(out.similarity,'VariableNames',{'CobraGurobi' 'MosekGurobi' 'MosekCobra'})
end