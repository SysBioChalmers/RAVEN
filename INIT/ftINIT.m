%This is the main ftINIT function
%Metabolomics is currently not supported, but will be when implemented into RAVEN.
function [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, taskReport, fullMipRes] = ftINIT(prepData, tissue, celltype, hpaData, transcrData, metabolomicsData, INITSteps, removeGenes, useScoresForTasks, paramsFT)
% ftINIT
%   Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks.
%   This is the newer version, which has updated handling of gene rules to
%   differentiate between isozymes and enzyme complexes.
%
%   prepData            The prepdata for the model.
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or transcrData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or transcrData.celltypes for this
%                       tissue (opt, default is to use the max values
%                       among all the cell types for the tissue.
%   hpaData             HPA data structure from parseHPA (opt if transcrData
%                       is supplied, default [])
%   transcrData         gene expression data structure (opt if hpaData is
%                       supplied, default []). Used to be called arrayData.
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not
%                       be unique, as there can be multiple cell types per
%                       tissue.
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%       threshold       a single value or a vector of gene expression 
%                       thresholds, above which genes are considered to be
%                       "expressed". default = 1(opt, by default, the mean expression
%                       levels of each gene across all tissues in transcrData
%                       will be used as the threshold values)
%       singleCells     binary value selecting whether to use the
%                       single-cell algorithm to identify expressed genes.
%                       If used, specify cell subpopulations in CELLTYPES
%                       (opt, default [])
%       plotResults     true if single cell probability distributions
%                       should be plotted (opt, default = false)
%   metabolomicsData    cell array with metabolite names that the model
%                       should produce (opt, default [])
%   rxnsToIgnorePatternStep1 Pattern describing which reactions to ignore in step 1 and 2:
%                       [b1,b2,b3,b4,b5,b6,b7,b8], bx is either 0 or 1, where 1 means that the group is excluded.
%                       b1 - Exchange rxns
%                       b2 - Import rxns without GPRs (from s into the cell)
%                       b3 - Simple transport reactions without GPRs (moves one metabolite between compartments)
%                       b4 - Advanced transport reactions without GPRs (moves metabolites between compartments, more complex function such as antiporter)
%                       b5 - Spontaneous reactions
%                       b6 - Reactions in the s compartment without GPRs
%                       b7 - Customly specified rxns (sent in when generating prepData)
%                       b8 - All rxns without GPRs
%   rxnsToIgnorePatternStep3 - same as above, but for step 3. Step 3 is only run if this differs from rxnsToIgnorePatternStep1
%   removeGenes         if true, low-abundance genes will be removed from
%                       grRules, unless they are the only gene associated 
%                       with a reaction, or a subunit of an enzyme complex
%                       (see "removeLowScoreGenes" function for details).
%                       If false, grRules will not be modified; however,
%                       genes that were associated only with removed 
%                       reactions will not be present in the final model.
%                       (opt, default true).
%   useScoresForTasks   true if the calculated reaction scored should be 
%                       used as weights when fitting to tasks (opt, default
%                       true)
%   milpSkipMets        Defines metabolites that will be removed from the model before running ftINIT (for example water)
%   allowExcretion      Allows excretion of metabolites in step 1 and 2.
%   printReport         true if a report should be printed to the screen
%                       (opt, default true)
%   params              parameter structure as used by getMILPParams. This
%                       is for the INIT algorithm. For the the MILP 
%                       problems solved to fit tasks, see paramsFT (opt,
%                       default [])
%   paramsFT            parameter structure as used by getMILPParams. This
%                       is for the fitTasks step. For the INIT algorithm,
%                       see params (opt, default [])
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
%   taskReport              structure with the results for each task
%   	id                  cell array with the id of the task
%       description         cell array with the description of the task
%       ok                  boolean array with true if the task was successful
%       essential           cell array with cell arrays of essential
%                           reactions for the task
%       gapfill             cell array of cell arrays of reactions included
%                           in the gap-filling for the task
%   fullMipRes              The solver results from the last MILP step run
%
%   This is the main function for automatic reconstruction of models based
%   on the INIT algorithm (PLoS Comput Biol. 2012;8(5):e1002518). Not all
%   settings are possible using this function, and you may want to call the
%   functions scoreComplexModel, runINIT and fitTasks individually instead.
%
%   NOTE: Exchange metabolites should normally not be removed from the model
%   when using this approach, since checkTasks/fitTasks rely on putting specific
%   constraints for each task. The INIT algorithm will remove exchange metabolites
%   if any are present. Use importModel(file,false) to import a model with
%   exchange metabolites remaining.
%
%   Usage: [model, metProduction, essentialRxnsForTasks, addedRxnsForTasks,...
%               deletedDeadEndRxns, deletedRxnsInINIT, taskReport] = ...
%               getINITModel2(refModel, tissue, celltype, hpaData, transcrData,...
%               metabolomicsData, removeGenes, taskFile, useScoresForTasks, ...
%               printReport, taskStructure, params, paramsFT);
%


if nargin < 5
    transcrData = [];
end
if nargin < 6
    metabolomicsData = [];
end
if nargin < 7
    INITSteps = getINITSteps([],'default');
end
if nargin < 8 || isempty(removeGenes)
    removeGenes = true;
end
if nargin < 9 || isempty(useScoresForTasks)
    useScoresForTasks = true;
end
if nargin < 10
    paramsFT = [];
end

%Handle detected mets:
%Previously, this was handled by giving a bonus for secreting those metabolites,
%but that doesn't work since the metabolite secretion and uptake can be lost when 
%we merge linearly dependent reactions.
%Instead, we need to figure out which reactions either produce or take up the mets.
%We then give a bonus if any of them carry flux.
%To simplify things, we focus on reactions that produce the metabolite (since there must be one such reaction). 
%It is still a bit complicated though. In this step, we focus on identifying
%producer reactions. We further reason that the direction doesn't matter - 
%we can force one of these reactions in any direction - if it becomes a consumer, it will
%automatically force another producer on as well (otherwise we'll have a net consumption).

if (~isempty(metabolomicsData))
    if length(unique(upper(metabolomicsData))) ~= length(metabolomicsData)
        dispEM('Metabolomics contains the same metabolite multiple times');
    end
    metData = false(numel(metabolomicsData), length(prepData.minModel.rxns)); %one row per metabolite that is a boolean vector
    for i=1:numel(metabolomicsData)
        %Get the matching mets
        metSel = ismember(upper(prepData.refModel.metNames),upper(metabolomicsData{i}));
        prodRxnsSel = any(prepData.refModel.S(metSel,:) > 0,1) | ... %direct producers
                     (any(prepData.refModel.S(metSel,:) < 0,1) & prepData.refModel.rev.'); %reversible reactions that are consumers
        %convert the production rxns from refModel to minModel
        prepData.groupIds
        [~,ia,ib] = intersect(prepData.minModel.rxns,prepData.refModel.rxns);
        grpIdsMerged = nan(length(prepData.minModel.rxns),1);
        grpIdsMerged(ia) = prepData.groupIds(ib);
        
        groupIdsPos = unique(prepData.groupIds(prodRxnsSel));%gets all group ids which includes a production rxn
        groupIdsPos = groupIdsPos(groupIdsPos ~= 0);%remove the 0 id, it means there is no group
        %the other option is that there is a direct match between the rxn id in minModel and refModel:
        posRxns = prepData.refModel.rxns(prodRxnsSel);
        directMatch = ismember(prepData.minModel.rxns, posRxns).';
        
        metData(i,:) = ismember(grpIdsMerged, groupIdsPos).' | directMatch;
    end
    metData = sparse(metData);
else
    metData = [];
end
%}

% Get rxn scores and adapt them to the minimized model
origRxnScores = scoreComplexModel(prepData.refModel,hpaData,transcrData,tissue,celltype);
origRxnScores(origRxnScores > -0.1 & origRxnScores <= 0) = -0.1;%we don't want reaction scores that are exactly 0 (or close), this causes problems in the milp
origRxnScores(origRxnScores < 0.1 & origRxnScores > 0) = 0.1;

rxnsTurnedOn = false(length(prepData.minModel.rxns),1);
fluxes = zeros(length(prepData.minModel.rxns),1);

rxnsToIgnoreLastStep = [1;1;1;1;1;1;1;1];

for initStep = 1:length(INITSteps)
    disp(['ftINIT: Running step ' initStep])
    stp = INITSteps{initStep};
    
    if any ((rxnsToIgnoreLastStep - stp.RxnsToIgnoreMask) < 0)
        dispEM('RxnsToIgnoreMask may not cover rxns not covered in previous steps, but the other way around is fine.');
    end
    rxnsToIgnoreLastStep = stp.RxnsToIgnoreMask;
    
    mm = prepData.minModel;
    
    if (~isempty(stp.MetsToIgnore))
        if (~isempty(milpSkipMets.simpleMets))
            %Here, we remove simple metabolites that will not really affect the milp but 
            %are very common in the S matrix. For example H2O, H+, etc.
            %It is also possible to leave compartments untouched, for example the i compartment in the mitochondria (for H+).
            metsToRem = ismember(mm.metNames,milpSkipMets.simpleMets.mets);
            compsToKeep = find(ismember(mm.comps, milpSkipMets.simpleMets.compsToKeep));
            metsToRem = metsToRem & ~ismember(mm.metComps, compsToKeep);
            mm.S(metsToRem,:) = 0;
        end    
    end

    %Set up the reaction scores and essential rxns
    rxnsToIgnore = getRxnsFromPattern(stp.RxnsToIgnoreMask, prepData);
    rxnScores = groupRxnScores(prepData.minModel, origRxnScores, prepData.refModel.rxns, prepData.groupIds, rxnsToIgnore);

    essentialRxns = prepData.essentialRxns;
    
    %Handle the results from previous steps ('ignore', 'exclude', 'essential')
    if strcmp(stp.HowToUsePrevResults, 'exclude')
        rxnScores(rxnsTurnedOn) = 0;
    elseif strcmp(stp.HowToUsePrevResults, 'essential')
        %Make all reversible reactions turned on in previous steps reversible
        %in the direction that they were previously carrying flux
        
        %first reverse the reactions that need to be reversed
        rev = mm.rev == 1;
        toRev = rxnsTurnedOn & rev & fluxes < 0;
        mm = reverseRxns(mm, mm.rxns(toRev));
        
        %Then make them irreversible
        mm.rev(toRev) = 0;
        mm.lb(toRev) = 0;

        essentialRxns = unique([prepData.essentialRxns;mm.rxns(rxnsTurnedOn)]);
    end

    
    mipGap = 1;
    first = true;
    success = false;
    for rn = 1:length(stp.MILPParams)
        params = stp.MILPParams{rn};
        if ~isfield(params, 'MIPGap')
            params.MIPGap = 0.0004;
        end
        
        if ~isfield(params, 'TimeLimit')
            params.TimeLimit = 5000;
        end
        
        if ~first 
            %There is sometimes a problem with that the objective function becomes close to zero,
            %which leads to that a small percentage of that (which is the MIPGap sent in) is very small
            %and the MILP hence takes a lot of time to finish. We also therefore use an absolute MIP gap, 
            %converted to a percentage using the last value of the objective function.
            params.MIPGap = min(max(params.MIPGap, stp.AbsMIPGaps{rn}/abs(lastObjVal)),1);
            params.seed = 1234;%use another seed, may work better

            if mipGap <= params.MIPGap
                success = true;
                break; %we're done - this will not happen the first time
            else
                disp(['MipGap too high, trying with a different run. MipGap = ' num2str(mipGap) ' New MipGap Limit = ' num2str(params.MIPGap)])
            end
        end
        
        first = false;
        
        %now run the MILP
        try
            %The prodweight for metabolomics is currently set to 5 - 0.5 was default in the old version, which I deemed very small?
            %There could be a need to specify this somewhere in the call at some point. 
            %This value has not been evaluated, but is assumed in the test cases - if changed, update the test case
            [deletedRxnsInINIT1, metProduction,fullMipRes,rxnsTurnedOn1,fluxes1] = ftINITInternalAlg(mm,rxnScores,metData,essentialRxns,5,stp.AllowMetSecr,stp.PosRevOff,params);
            mipGap = fullMipRes.mipgap;
            lastObjVal = fullMipRes.obj;
        catch e
            mipGap = Inf;
            lastObjVal = Inf; %we need to set something here, Inf leads to that this doesn't come into play
        end
        
        success = mipGap <= params.MIPGap;
    end
    
    if ~success
        dispEM(['Failed to find good enough solution within the time frame. MIPGap: ' num2str(mipGap)]);
    end
    
    %save the reactions turned on and their fluxes for the next step
    rxnsTurnedOn = rxnsTurnedOn | rxnsTurnedOn1.';
    %The fluxes are a bit tricky - what if they change direction between the steps?
    %The fluxes are used to determine the direction in which reactions are forced on 
    %(to simplify the problem it is good if they are unidirectional).
    %We use the following strategy:
    %1. Use the fluxes from the most recent step.
    %2. If any flux is very low there (i.e. basically zero), use the flux from the previous steps
    %This could in theory cause problems, but seems to work well practically
    fluxesOld = fluxes;
    fluxes = fluxes1;
    fluxes(abs(fluxes1) < 10^-7) = fluxesOld(abs(fluxes1) < 10^-7);
end


%get the essential rxns
essential = ismember(prepData.minModel.rxns,prepData.essentialRxns);
%So, we only add reactions where the linearly merged scores are zero for all linearly dependent reactions 
% (this cannot happen by chance, taken care of in the function groupRxnScores)
rxnsToIgn = rxnScores == 0; 
deletedRxnsInINITSel = ~(rxnsTurnedOn | rxnsToIgn | essential);
deletedRxnsInINIT = prepData.minModel.rxns(deletedRxnsInINITSel);

%Here we need to figure out which original reactions (before the linear merge) 
%that were removed. These are all reactions with the same group ids as the removed reactions
groupIdsRemoved = prepData.groupIds(ismember(prepData.refModel.rxns, deletedRxnsInINIT)); %can improve this slightly, use sel above
groupIdsRemoved = groupIdsRemoved(groupIdsRemoved ~= 0);%zero means that the reaction was not grouped, all with zeros are not a group!
rxnsToRem = union(prepData.refModel.rxns(ismember(prepData.groupIds,groupIdsRemoved)), deletedRxnsInINIT);%make a union here to include the ungrouped (unmerged) as well

initModel = removeReactions(prepData.refModel,rxnsToRem,false,true);

% remove metabolites separately to avoid removing those needed for tasks
unusedMets = initModel.mets(all(initModel.S == 0,2));
initModel = removeMets(initModel, setdiff(unusedMets, prepData.essentialMetsForTasks));

%if printReport == true
%    printScores(initModel,'INIT model statistics',hpaData,transcrData,tissue,celltype);
%    printScores(removeReactions(cModel,setdiff(cModel.rxns,rxnsToRem),true,true),'Reactions deleted by INIT',hpaData,transcrData,tissue,celltype);
%end

%The full model has exchange reactions in it. fitTasks calls on fillGaps,
%which automatically removes exchange metabolites (because it assumes that
%the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are deleted from the
%reference model here.
initModel.id = 'INITModel';

%If gaps in the model should be filled using a task list
if ~isempty(prepData.taskStruct)
    %Remove exchange reactions and reactions already included in the INIT
    %model
    %We changed strategy and instead include all rxns except the exchange rxns in the ref model
    %But we do keep the exchange rxns that are essential.
    %Let's test to remove all, that should work
    
    %At this stage the model is fully connected and most of the genes with
    %good scores should have been included. The final gap-filling should
    %take the scores of the genes into account, so that "rather bad"
    %reactions are preferred to "very bad" reactions. However, reactions
    %with positive scores will be included even if they are not connected
    %in the current formulation. Therefore, such reactions will have to be
    %assigned a small negative score instead.
    exchRxns = getExchangeRxns(prepData.refModel);
    refModelNoExc = removeReactions(prepData.refModelWithBM,exchRxns,true,true);
    exchRxns = getExchangeRxns(initModel);
    initModelNoExc = removeReactions(addBoundaryMets(initModel),exchRxns,true,true);
    
    if useScoresForTasks == true
        %map the rxn scores to the model without exchange rxns
        [~,ia,ib] = intersect(refModelNoExc.rxns,prepData.refModel.rxns);
        rxnScores2nd = NaN(length(refModelNoExc.rxns),1);
        rxnScores2nd(ia) = origRxnScores(ib);
        %all(rxnScores2nd == refRxnScores);%should be the same, ok!
        %temp for testing
        global g_initModel g_refModelNoExc g_refRxnScores g_taskStructure g_paramsFT; 
        g_initModel = initModel;
        g_refModelNoExc = refModelNoExc;
        g_refRxnScores = rxnScores2nd;
        g_taskStructure = prepData.taskStruct;
        g_paramsFT = paramsFT; 
        [outModel,addedRxnMat] = fitTasksOpt(initModelNoExc,refModelNoExc,[],true,min(rxnScores2nd,-0.1),prepData.taskStruct,paramsFT);
    else
        [outModel,addedRxnMat] = fitTasksOpt(initModelNoExc,refModelNoExc,[],true,[],prepData.taskStruct,paramsFT);
    end
    %if printReport == true
    %    printScores(outModel,'Functional model statistics',hpaData,transcrData,tissue,celltype);
    %    printScores(removeReactions(outModel,intersect(outModel.rxns,initModel.rxns),true,true),'Reactions added to perform the tasks',hpaData,transcrData,tissue,celltype);
    %end
    
    addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));
else
    outModel = initModel;
    addedRxnMat = [];
    addedRxnsForTasks = {};
end

% The model can now perform all the tasks defined in the task list.
model = outModel;


% At this stage the model will contain some exchange reactions but probably
% not all (and maybe zero). This can be inconvenient, so all exchange
% reactions from the reference model are added, except for those which
% involve metabolites that are not in the model.



%TODO: NOT SURE ABOUT THE TASK REPORT THINGS BELOW, INVESTIGATE WHAT THIS IS USED FOR
% The task analysis should probably be done separately, for implementation in RAVEN, remove this part.
% Add information about essential reactions and reactions included for
% gap-filling and return a taskReport
taskReport = prepData.taskReport;
if ~isempty(prepData.taskStruct)
    I = find(taskReport.ok); %Ignore failed tasks
    for i = 1:numel(I)
        taskReport.gapfill{I(i),1} = refModelNoExc.rxns(addedRxnMat(:,i));
    end
else
    taskReport = [];
end


%Start from the original model, and just remove the reactions that are no longer there (and keep exchange rxns). The model we got out
%from the problem is not complete, it doesn't have GRPs etc.
%The logic below is a bit complicated. We identify the reactions that should be removed from the full model as 
%reactions that have been removed in the init model except the ones that were added back. In addition, we make 
%sure that no exchange rxns are removed - they can be removed in the init model if they were linearly merged with other
%reactions that were decided to be removed from the model. We want to keep all exchange rxns to make sure the tasks can
%be performed also without manipulating the b vector in the model (which is what is done in the gap-filling).
exchRxns = getExchangeRxns(prepData.refModel);
deletedRxnsInINIT = setdiff(prepData.refModel.rxns,union(union(initModel.rxns, addedRxnsForTasks), exchRxns));
outModel = removeReactions(prepData.refModel, deletedRxnsInINIT, true); %we skip removing the genes for now, I'm not sure it is desirable

% If requested, attempt to remove negative-score genes from the model, 
% depending on their role (isozyme or complex subunit) in each grRule.
% See the "removeLowScoreGenes" function more more details, and to adjust
% any default parameters therein.
if ( removeGenes )
    [~, geneScores] = scoreComplexModel(outModel,hpaData,transcrData,tissue,celltype);
    outModel = removeLowScoreGenes(outModel,geneScores);
end


model = outModel;

end

%This is for printing a summary of a model
function [rxnS, geneS] = printScores(model,name,hpaData,transcrData,tissue,celltype)
[a, b] = scoreComplexModel(model,hpaData,transcrData,tissue,celltype);
rxnS = mean(a);
geneS = mean(b,'omitnan');
fprintf([name ':\n']);
fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end

function rxnsToIgnore = getRxnsFromPattern(rxnsToIgnorePattern, prepData)
    rxnsToIgnore = false(length(prepData.toIgnoreExch),1);
    if rxnsToIgnorePattern(1) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreExch; end;
    if rxnsToIgnorePattern(2) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreImportRxns; end;
    if rxnsToIgnorePattern(3) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreSimpleTransp; end;
    if rxnsToIgnorePattern(4) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreAdvTransp; end;
    if rxnsToIgnorePattern(5) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreSpont; end;
    if rxnsToIgnorePattern(6) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreS; end;
    if rxnsToIgnorePattern(7) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreCustomRxns; end;
    if rxnsToIgnorePattern(8) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreAllWithoutGPRs; end;
end

