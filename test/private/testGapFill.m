function out=testGapFill(param)
%if (nargin<1 || isempty(param))
%	param.solver='gurobi';
%	param.biomassReactionExcludeRatio=.4;
%	param.useConstraints=true;
%end
def.solver='gurobi';
def.biomassReactionExcludeRatio=.4;
def.useConstraints=true;

param=structUpdate(def,param);

modLoad('~all RAVENgit',[]);
modelPath='/Users/hedani/Documents/GitRepos/yeast-metabolic-network-7.6/';
model=importModel([modelPath 'yeast_7.6_cobra.xml'],true,true);
model=setParam(model,'obj',{'r_4041'},1); % biomass
model=setParam(model,'lb',{'r_4041'},.1); % biomass

%model.lb(end)=.1;

refModel=model;

SF=full(model.S);
biomassMetaboliteInd=find(abs(SF(:,end))>0);
biomassReactionInd=arrayfun(@(x) find(abs(SF(x,1:end-1))>0),biomassMetaboliteInd,'UniformOutput',false);
biomassReactionInd=unique([biomassReactionInd{:}]);
%biomassReactionInd=setdiff(biomassReactionInd,{4041});
biomassReactionNames=model.rxns(biomassReactionInd);

sel=randsample([1:length(biomassReactionNames)],floor(param.biomassReactionExcludeRatio*length(biomassReactionNames)));

biomassReactionSubsetInd=biomassReactionInd(sel);
biomassReactionSubsetNames=biomassReactionNames(sel);
%biomassReactionSubsetInd=[2064];
%biomassReactionSubsetNames=biomassReactionNames([2064]);

%model=removeReactions(model,randi([1,length(model.rxns)],1,200));

model=removeReactions(model,biomassReactionSubsetInd); %biomass metabolite heavy reaction set
model.id='yeast_reduced';
solveLP(model)
%disp('Using mosek')
%setRavenSolver('mosek')
%[newConnected cannotConnect addedRxns newModel exitFlag]=fillGaps(model,refModel)
%mosekRes={newConnected cannotConnect addedRxns newModel exitFlag};

setRavenSolver(param.solver)
[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(model,refModel,false,param.useConstraints);
gurobiRes={newConnected cannotConnect addedRxns newModel exitFlag};


moreRes={biomassReactionNames,biomassReactionInd,biomassReactionSubsetNames,biomassReactionSubsetInd};
%disp(addedRxns)
%disp(biomassReactionSubsetNames)


foundI=cellfun(@(x) find(strcmp(biomassReactionSubsetNames,x)),addedRxns,'UniformOutput',false);
foundI=[foundI{:}];
nonfoundI=setdiff([1:length(biomassReactionSubsetNames)],foundI);
disp('Added biomass metabolite producing/consuming reactions:')
disp(sort({biomassReactionSubsetNames{foundI}}))

disp('Biomass metabolite producing/consuming reactions not added:')
disp(sort({biomassReactionSubsetNames{nonfoundI}}))

%disp('Biomass metabolite producing/consuming reactions total removed:')
%disp(sort(biomassReactionSubsetNames)')
%cellfun(@(x, y) disp([x '      ' y]),biomassReactionSubsetNames{foundI},'UniformOutput',false,'ErrorHandler',@(S,varargin) disp(varargin))
%cellfun(@(x) disp(x),biomassReactionSubsetNames{nonfoundI},'UniformOutput',false)
solveLP(newModel)
out={gurobiRes,moreRes,model};

% check if all biomass stuff added



end

function s_merged=structUpdate(s_old,s_new)
s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
names = [fieldnames(s_merged); fieldnames(s_new)];
s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end
