function [model, metProduction, essentialRxnsForTasks, addedRxnsForTasks, deletedDeadEndRxns, deletedRxnsInINIT, taskReport]=getINITModel(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT)
% getINITModel
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%
%   refModel            a model structure. The model should be in the
%                       closed form (no exchange reactions open). Import
%                       using import(filename,false). If the model is not
%                       loaded using importModel, it might be that there
%                       is no "unconstrained" field. In that case,
%                       manually add the field like:
%                       model.unconstrained=false(numel(model.mets),1);
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this
%                       tissue (opt, default is to use the best values
%                       among all the cell types for the tissue. Use [] if
%                       you want to supply more arguments)
%   hpaData             HPA data structure from parseHPA (opt if arrayData is
%                       supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not be
%                       unique, as there can be multiple cell types per tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%   metabolomicsData    cell array with metabolite names that the model
%                       should produce (opt, default [])
%   taskFile            a task list in Excel format. See parseTaskList for
%                       details (opt, default [])
%   useScoresForTasks   true if the calculated reaction scored should be used as
%                       weights in the fitting to tasks (opt, default true)
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   taskStructure       task structure as from parseTaskList. Can be used
%                       as an alternative way to define tasks when Excel
%                       sheets are not suitable. Overrides taskFile (opt,
%                       default [])
%   params              parameter structure as used by getMILPParams. This is
%                       for the INIT algorithm. For the the MILP problems
%                       solved to fit tasks, see paramsFT (opt, default [])
%   paramsFT            parameter structure as used by getMILPParams. This is
%                       for the fitTasks step. For the INIT algorithm, see
%                       params (opt, default [])
%
%   model                   the resulting model structure
%   metProduction           array that indicates which of the
%                           metabolites in metabolomicsData that could be
%                           produced. Note that this is before the
%                           gap-filling process to enable defined tasks. To
%                           see which metabolites that can be produced in
%                           the final model, use canProduce.
%                           -2: metabolite name not found in model
%                           -1: metabolite found, but it could not be produced
%                           1: metabolite could be produced
%   essentialRxnsForTasks   cell array of the reactions which were
%                           essential to perform the tasks
%   addedRxnsForTasks       cell array of the reactions which were added in
%                           order to perform the tasks
%   deletedDeadEndRxns      cell array of reactions deleted because they
%                           could not carry flux (INIT requires a
%                           functional input model)
%   deletedRxnsInINIT       cell array of the reactions which were deleted by
%                           the INIT algorithm
%   taskReport              structure with the results for each task
%   	id                  cell array with the id of the task
%       description         cell array with the description of the task
%       ok                  boolean array with true if the task was successful
%       essential           cell array with cell arrays of essential
%                           reactions for the task
%       gapfill             cell array of cell arrays of reactions included
%                           in the gap-filling for the task
%
%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreModel, runINIT and fitTasks individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasks rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [model, metProduction, essentialRxnsForTasks, addedRxnsForTasks,...
%               deletedDeadEndRxns, deletedRxnsInINIT, taskReport]=...
%               getINITModel(refModel, tissue, celltype, hpaData, arrayData,...
%               metabolomicsData, taskFile, useScoresForTasks, printReport,...
%               taskStructure, params, paramsFT)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<3
    celltype=[];
end
if nargin<4
    hpaData=[];
end
if nargin<5
    arrayData=[];
end
if nargin<6
    metabolomicsData=[];
end
if nargin<7
    taskFile=[];
end
if nargin<8
    useScoresForTasks=true;
end
if nargin<9
    printReport=true;
end
if nargin<10
    taskStructure=[];
end
if nargin<11
    params=[];
end
if nargin<12
    paramsFT=[];
end

%Check that the model is in the closed form
if ~isfield(refModel,'unconstrained')
    EM='Exchange metabolites should normally not be removed from the model when using getINITModel. Use importModel(file,false) to import a model with exchange metabolites remaining (see the documentation for details)';
    dispEM(EM);
end

%Create the task structure if not supplied
if any(taskFile) && isempty(taskStructure)
	taskStructure=parseTaskList(taskFile);
end

if printReport==true
    if any(celltype)
        fprintf(['***Generating model for: ' tissue ' - ' celltype '\n']);
    else
        fprintf(['***Generating model for: ' tissue '\n']);
    end
    if ~isempty(hpaData)
        fprintf('-Using HPA data\n');
    end
    if ~isempty(arrayData)
        fprintf('-Using array data\n');
    end
    if ~isempty(metabolomicsData)
        fprintf('-Using metabolomics data\n');
    end
    if any(taskFile)
        fprintf('-Using metabolic tasks\n');
    end
    fprintf('\n');

    printScores(refModel,'Reference model statistics',hpaData,arrayData,tissue,celltype);
end

%Remove dead-end reactions to speed up the optimization and to
%differentiate between reactions removed by INIT and those that are
%dead-end
[~, deletedDeadEndRxns]=simplifyModel(refModel,true,false,true,true,true);
cModel=removeReactions(refModel,deletedDeadEndRxns,false,true);

%Store the connected model like this to keep track of stuff
if printReport==true
    printScores(cModel,'Pruned model statistics',hpaData,arrayData,tissue,celltype);
end

%If tasks have been defined, then go through them and get essential
%reactions
if ~isempty(taskStructure)
    [taskReport, essentialRxnMat]=checkTasks(cModel,[],printReport,true,true,taskStructure);

    essentialRxnsForTasks=cModel.rxns(any(essentialRxnMat,2));

    %Remove tasks that cannot be performed
    taskStructure(taskReport.ok==false)=[];
    if printReport==true
        printScores(removeReactions(cModel,setdiff(cModel.rxns,essentialRxnsForTasks),true,true),'Reactions essential for tasks',hpaData,arrayData,tissue,celltype);
    end
else
    essentialRxnsForTasks={};
end

%Score the connected model
[rxnScores, geneScores]=scoreModel(cModel,hpaData,arrayData,tissue,celltype);

%Run the INIT algorithm. The exchange reactions that are used in the final
%reactions will be open, which doesn't fit with the last step. Therefore
%delete reactions from the original model instead of taking the output.
%The default implementation does not constrain reversible reactions to only
%carry flux in one direction.
%Runs without the constraints on reversibility and with all output allowed.
%This is to reduce the complexity of the problem.
[~, deletedRxnsInINIT, metProduction]=runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,true,false,params);
initModel=removeReactions(cModel,deletedRxnsInINIT,true,true);
if printReport==true
    printScores(initModel,'INIT model statistics',hpaData,arrayData,tissue,celltype);
    printScores(removeReactions(cModel,setdiff(cModel.rxns,deletedRxnsInINIT),true,true),'Reactions deleted by INIT',hpaData,arrayData,tissue,celltype);
end

%The full model has exchange reactions in it. fitTasks calls on fillGaps,
%which automatically removes exchange metabolites (because it assumes that
%the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are delete from the
%reference model here.
initModel.id='INITModel';

%If gaps in the model should be filled using a task list
if ~isempty(taskStructure)
    %Remove exchange reactions and reactions already included in the INIT
    %model
    refModelNoExc=removeReactions(refModel,union(initModel.rxns,getExchangeRxns(refModel)),true,true);

    %At this stage the model is fully connected and most of the genes with
    %good scores should have been included. The final gap-filling should
    %take the scores of the genes into account, so that "rather bad"
    %reactions are preferred to "very bad" reactions. However, reactions
    %with positive scores will be included even if they are not connected
    %in the current formulation. Therefore, such reactions will have to be
    %assigned a small negative score instead.
    if useScoresForTasks==true
        refRxnScores=scoreModel(refModelNoExc,hpaData,arrayData,tissue,celltype);
        [outModel, addedRxnMat]=fitTasks(initModel,refModelNoExc,[],true,min(refRxnScores,-0.1),taskStructure,paramsFT);
    else
        [outModel, addedRxnMat]=fitTasks(initModel,refModelNoExc,[],true,[],taskStructure,paramsFT);
    end
    if printReport==true
        printScores(outModel,'Functional model statistics',hpaData,arrayData,tissue,celltype);
        printScores(removeReactions(outModel,intersect(outModel.rxns,initModel.rxns),true,true),'Reactions added to perform the tasks',hpaData,arrayData,tissue,celltype);
    end

    addedRxnsForTasks=refModelNoExc.rxns(any(addedRxnMat,2));
else
    outModel=initModel;
    addedRxnMat=[];
    addedRxnsForTasks={};
end

%The model can now perform all the tasks defined in the task list. The
%algorithm cannot deal with gene-complexes at the moment. It is therefore
%ok to remove bad genes from a reaction (as long as at least one gene is
%kept)
model=outModel;

[~, I]=ismember(model.genes,cModel.genes); %All should be found
%This is a little weird way to make sure that only one bad gene is included
%if there are no good ones (since all -Inf==max(-Inf))
geneScores(isinf(geneScores))=-1000+rand(sum(isinf(geneScores)),1);

model.grRules(:)={''};
for i=1:numel(model.rxns)
    ids=find(model.rxnGeneMat(i,:));
    if numel(ids)>1
       scores=geneScores(I(ids));
       %Only keep the positive ones if possible
       model.rxnGeneMat(i,ids(~(scores>0 | scores==max(scores))))=0;
    end
    %Rewrite the grRules to be only OR
    if isfield(model,'grRules')
       J=find(model.rxnGeneMat(i,:));
       for j=1:numel(J)
           model.grRules{i}=[model.grRules{i} '(' model.genes{J(j)} ')'];
           if j<numel(J)
               model.grRules{i}=[model.grRules{i} ' or '];
           end
       end
    end
end

%Find all genes that are not used and delete them
I=sum(model.rxnGeneMat)==0;
model.genes(I)=[];
model.rxnGeneMat(:,I)=[];
if isfield(model,'geneShortNames')
    model.geneShortNames(I)=[];
end
if isfield(model,'geneMiriams')
    model.geneMiriams(I)=[];
end
if isfield(model,'geneFrom')
    model.geneFrom(I)=[];
end
if isfield(model,'geneComps')
    model.geneComps(I)=[];
end

%At this stage the model will contain some exchange reactions but probably not all
%(and maybe zero). This can be inconvenient, so all exchange reactions from the
%reference model are added, except for those which involve metabolites that
%are not in the model.

%First delete and included exchange reactions in order to prevent the order
%from changing
model=removeReactions(model,getExchangeRxns(model));

%Create a model with only the exchange reactions in refModel
excModel=removeReactions(refModel,setdiff(refModel.rxns,getExchangeRxns(refModel)),true,true);

%Find the metabolites there which are not exchange metabolites and which do
%not exist in the output model
I=~ismember(excModel.mets,model.mets) & excModel.unconstrained==0;

%Then find those reactions and delete them
[~, J]=find(excModel.S(I,:));
excModel=removeReactions(excModel,J,true,true);

%Merge with the output model
model=mergeModels({model;excModel});
model.id='INITModel';
model.description=['Automatically generated model for ' tissue];
if any(celltype)
    model.description=[model.description ' - ' celltype];
end

if printReport==true
    printScores(model,'Final model statistics',hpaData,arrayData,tissue,celltype);
end

%Add information about essential reactions and reactions included for
%gap-filling and return a taskReport
if ~isempty(taskStructure)
    I=find(taskReport.ok); %Ignore failed tasks
    for i=1:numel(I)
        taskReport.essential{I(i),1}=cModel.rxns(essentialRxnMat(:,I(i)));
        taskReport.gapfill{I(i),1}=refModelNoExc.rxns(addedRxnMat(:,i));
    end
else
    taskReport=[];
end
end

%This is for printing a summary of a model
function [rxnS, geneS]=printScores(model,name,hpaData,arrayData,tissue,celltype)
    [a, b]=scoreModel(model,hpaData,arrayData,tissue,celltype);
    rxnS=mean(a);
    geneS=mean(b(~isinf(b)));
    fprintf([name ':\n']);
    fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
    fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
    fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
    fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end
