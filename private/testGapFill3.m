function out=testGapFill2(param)
	if nargin<1
		param.solver='gurobi';
		param.useConstraints=true;
	end

modelPath='/Users/hedani/Documents/GitRepos/yeast-metabolic-network-7.6/';
model=importModel([modelPath 'yeast_7.6_cobra.xml'],true,true);

IndexC = regexpi(model.rxnNames, 'glucose');
Index = find(not(cellfun('isempty', IndexC)));

setRavenSolver(param.solver)

nonEssential={};
essential={};
models={};

for i=1:length(Index)
	model2=removeReactions(model,Index(i));
	model2=setParam(model2,'obj',{'r_4041'},1); % biomass
	model2=setParam(model2,'lb',{'r_4041'},.1);
	model2.id='yeast_reduced';
	[newConnected cannotConnect addedRxns newModel exitFlag]=fillGaps(model2,model,false,true)
	if (length(addedRxns)==0)
		nonEssential{end+1}=model.rxns{Index(i)};
	else
		essential{end+1}=model.rxns{Index(i)};
	end
	models{end+1}=newModel;
	disp('Test model feasibility')
	solveLP(newModel)
end

model2=removeReactions(model,Index,true);
model2=setParam(model2,'obj',{'r_4041'},1); % biomass
model2=setParam(model2,'lb',{'r_4041'},.1);
model2.id='yeast_reduced';
[newConnected cannotConnect addedRxns newModel exitFlag]=fillGaps(model2,model,false,true)

%solveLP(newModel)

%newModel=removeReactions(newModel,{'r_1084','r_0466'});

%solveLP(newModel)

disp('Essential rxns:')
cellfun(@(rxn) printModel(model,rxn),essential)
disp('Non essential rxns:')
cellfun(@(rxn) printModel(model,rxn),nonEssential)

out=newModel;

end

% >> model.rxns(Index)

% ans = 

%     'r_0466'
%     'r_0467'
%     'r_0534'
%     'r_1068'
%     'r_1070'
%     'r_1071'
%     'r_1084'
%     'r_1166'
%     'r_1714'
%     'r_1805'