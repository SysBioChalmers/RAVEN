function prepData = prepINITModel(origRefModel, taskStruct, spontRxnNames, convertGenes, customRxnsToIgnore, extComp, skipScaling)
% prepINITModel
%
% The purpose of this function is to run time-consuming calculation steps that are not
% dependent on the RNA-Seq data.
%
% origRefModel      The model to use. Expected to be something such as Human-GEM, 
%                   Mouse-GEM, etc.
% taskStruct        The essential tasks. Can be loaded with for example
%                   taskStruct = parseTaskList('../data/metabolicTasks_Essential.txt');
% spontRxnNames     The spontaneous rxns. (optional, default {})
% convertGenes      If true the genes are converted to gene names (from 
%                   ENSEMBL) (optional, default false)
% customRxnsToIgnore These reactions can be ignored in the ignore mask 
%                   (specifying b7=1) (optional, default = {})
% extComp           Name of the external compartment, typically 's' or 'e'. This
%                   is used for identifying exch and import rxns (optional, default = 'e')
% prepData          The resulting prepData structure which is used as input to ftINIT
% skipScaling       If true the scaling step is not run on the minimal model. The 
%                   scaling is there to remove large differences between the 
%                   stoichiometric coefficients within a reaction, since such
%                   differences creates numerical issues in the ftINIT algorithm 
%                   due to limitations in solver resolution. However,
%                   the current scaling step is also risky and may lead to that
%                   some reactions cannot carry flux, since it changes the stoichiometry
%                   of the reactions with large differences. If you experience problems
%                   where the solution is infeasible, it may be worth trying to turn off
%                   the scaling. Note that it is only the minModel that is scaled,
%                   the scaling will not be present in the final model.
%                   Default: (optional, default = false)
% 
% Usage: prepData = prepINITModel(origRefModel, taskStruct, spontRxnNames, convertGenes, customRxnsToIgnore, extComp)


if nargin < 3
    spontRxnNames = {}; 
end

if nargin < 4
    convertGenes = false; 
end

if nargin < 5
    customRxnsToIgnore = {}; 
end

if nargin < 6
    extComp = 'e';
end

if nargin < 7
    skipScaling = false;
end
disp('Step 1: Gene rules')
[origRefModel.grRules, origRefModel.rxnGeneMat] = standardizeGrRules(origRefModel, true);

if convertGenes %For mouse we might want to translate in the opposite direction - this has to be done before calling this function in that case.
    [origRefModel.grRules, origRefModel.genes, origRefModel.rxnGeneMat] = translateGrRules(origRefModel.grRules, 'Name');
end


% Get list of all dead-end/inacessible/constrained reactions in the model.
% We don't merge linear dependent reactions here, that will happen later.
% We also don't remove reversibility here, we don't want that in the final
% models produced. Therefore, this is done as a separate step.
% This takes quite some time to run
disp('Step 2: First simplification')
[~,deletedDeadEndRxns] = simplifyModel(origRefModel,true,false,true,true,true);

% get reduced model
cModel = removeReactions(origRefModel,deletedDeadEndRxns,false,true);


%Then, run the tasks to find out which reactions are essential. These rxns
%will be forced to carry flux in the optimization later

% determine reactions essential for tasks - takes ~10 min
% need to add boundary mets first - but only to a temporary model, we don't
% want those for further processing
disp('Step 3: Check tasks (~10 min)')
if ~isempty(taskStruct)
    bModel = closeModel(cModel);
    [taskReport, essentialRxnMat, ~, essentialFluxes] = checkTasks(bModel,[],true,false,true,taskStruct);

    %extract the essential rxns:
    sel = sum(essentialRxnMat,2) > 0;
    selInd = find(sel);
    essentialRxns = bModel.rxns(selInd);%382 rxns for human-GEM

    %Find metabolites present in taskStruct. We want to avoid removing
    %these metabolites from the final model (even though they may no longer
    %participate in any reacitons) so that the final model is still able to
    %complete all of the tasks without any errors.
    taskMets = union(vertcat(taskStruct.inputs),vertcat(taskStruct.outputs));
    taskMets = union(taskMets, parseRxnEqu(vertcat(taskStruct.equations)));
    modelMets = strcat(cModel.metNames,'[',cModel.comps(cModel.metComps),']');
    [inModel,metInd] = ismember(taskMets,modelMets);
    essentialMetsForTasks = cModel.mets(metInd(inModel));%118 mets

    %Remove tasks that cannot be performed
    taskStruct(taskReport.ok == false) = [];

else
    essentialRxns = {};
    essentialMetsForTasks = {};
    taskReport = {};
end


% remove metabolites separately to avoid removing those needed for tasks
unusedMets = cModel.mets(all(cModel.S == 0,2));%1248 in Human-GEM
cModel = removeMets(cModel, setdiff(unusedMets,essentialMetsForTasks));

if ~isempty(taskStruct)

    %Here, we decide on a direction in which the essential reactions should be forced.
    essentialRevDir = false(length(essentialRxns),1);
    pp = zeros(length(essentialRxns),1);
    nn = zeros(length(essentialRxns),1);
    for i = 1:length(selInd)
        pos = sum(essentialFluxes(selInd(i),essentialRxnMat(selInd(i),:)) > 0);
        neg = sum(essentialFluxes(selInd(i),essentialRxnMat(selInd(i),:)) < 0);
        essentialRevDir(i) = pos < neg;
        pp(i) = pos;
        nn(i) = neg;
    end

    %sum(pp>0 & nn>0) %0, so it is not a big problem that different tasks share essential 
    %reactions but in different directions, at least no such cases exist now, so we ignore this for now.

    %Now create a minimal model for the MILP - we want to make this as small as we can, but we
    %don't want to use this later when we create the final models - just for the first MILP.
    %We don't want boundary metabolites in here either.

    %We do a trick with the essential rxns here:
    %What we want is to get rid of the reversibility of the reactions, because reversible reactions require
    %an integer each in the MILP, and we can also get the case that it may be cheaper to force flux in the wrong direction, which 
    %will probably not give the gap-filling of reactions that we want.
    %So, we simply change the reversible reactions to irreversible, and flip their direction if the direction
    %to force flux in is negative. We don't need to do anything else, runINIT will then handle all essential rxns
    %as irreversible (and will not add any integers etc. for them).

    minModel1 = cModel;
    %for the essential rxns where we should force flux in the positive direction, just make them irreversible
    %for the negative direction, we do the same, but also flip the direction of the reaction

    %first flip the reactions with negative direction
    %constructEquations(minModel1,minModel1.rxns(selInd(essentialRevDir))) %looks good
    minModel1 = reverseRxns(minModel1, minModel1.rxns(selInd(essentialRevDir)));
    %constructEquations(minModel2,minModel2.rxns(selInd(essentialRevDir))) %looks good

    %Then make them irreversible
    minModel1.rev(selInd) = 0;
    minModel1.lb(selInd) = 0;
else
    minModel1 = cModel;
end
%Now, remove reversibility on any reactions that can only be run in one direction - that simplifies the problem
%This takes a long time to run (at least an hour or so).
%This is an important step, we go from 5533 to 3734 reversible rxns, which reduces the size of the MILP,
%since irreversible rxns are easier to handle.
disp('Step 4: Second simplification (~1 hour)')
minModel2 = simplifyModel(minModel1,false,false,false,false,false,false,true);

disp('Step 5: Final work')

%Now, merge all linear dependent rxns - this reduces the size of the model a lot.
%Size goes from 11888 rxns to 7922 rxns for human-GEM, so this is an important step
[minModel3,origRxnIds,groupIds]=mergeLinear(minModel2, {});

%we now need to figure out what happened to the essential rxns in the merge step above.
%Note that it would be ideal to run checkTasks on the result from mergeLinear, but that doesn't work.
%The exchange rxns have in many cases been merged with other reactions, making it difficult to run
%the tasks. So, instead we run checkTasks before, and figure out which merged reactions the 
%essential reactions belong to now

%It may be possible to create a minModel2b, run mergeLinear on that (which is quick), and 
%then run checkTasks on the minimized model. Would be good to check if this gives the 
%same result, I'm not sure. The code would get less complicated.


if ~isempty(taskStruct)

    %find all potential rxns
    rxnIndOrig = ismember(origRxnIds, essentialRxns);
    ids = unique(groupIds(rxnIndOrig));
    ids = ids(ids > 0);%0 means it has not been merged
    rxnCandidates = unique([essentialRxns;origRxnIds(ismember(groupIds, ids))]);%388, essentialRxns is 382
    newEssentialRxns = minModel3.rxns(ismember(minModel3.rxns,rxnCandidates));%196 in Human-GEM
    %setdiff(newEssentialRxns, essentialRxns)%only 2, that is very few
    %setdiff(essentialRxns, newEssentialRxns)%188, that is very many
    %It is probably the case that many essential reactions are linearly dependent on one another, which makes sense.
else
    newEssentialRxns = {};
end
%What is left now is to remove some rxns from the problem, such as exchange rxns and spontaneous rxns.
%We set the rxn scores of these to exactly 0, and make sure no other scores are zero.
%runINIT will then handle this, we do >0 and <0 to get positive and negative scores, the rest 
%are not included in the problem.
%But, this has to happen in the scoring - we identify the reactions in the origrxns here, and then fix this
%in the scoring for each sample (the function groupRxnScores handles this and is run for each sample).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the reactions to ignore in the minModel2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Identify reactions that should be ignored by tINIT, i.e. tINIT will not remove these regardless of score


[~,exchRxnInd] = getExchangeRxns(minModel2);
toIgnoreExch = false(numel(minModel2.rxns),1);
toIgnoreExch(exchRxnInd) = true;
toIgnoreImportRxns = false(numel(minModel2.rxns),1);
toIgnoreSimpleTransp = false(numel(minModel2.rxns),1);
toIgnoreAdvTransp = false(numel(minModel2.rxns),1);
toIgnoreS = false(numel(minModel2.rxns),1);



%sum(toIgnore) %1500 for Human-GEM
%find simple transport reactions with no GPRs:
%We really only skip the reactions with two metabolites
%where they are the same but different compartments - only skip the transport into the cell for now
numMets = sum(minModel2.S ~= 0, 1);
scomp = find(strcmp(minModel2.comps,extComp));
for i = 1:length(minModel2.rxns)
    %check that it has no GPR and that the rxn has two metabolites
    if strcmp(minModel2.grRules(i),'') && (numMets(i) == 2)
      metsComps = minModel2.metComps(minModel2.S(:,i) ~= 0);
      metNames = minModel2.metNames(minModel2.S(:,i) ~= 0);
      if metsComps(1) ~= metsComps(2) && (metsComps(1) == scomp || metsComps(2) == scomp) %different compartments, one is the extComp
          if strcmp(metNames{1}, metNames{2}) %the same metabolite, it is a transport reaction
              toIgnoreImportRxns(i) = true;
          end
      elseif metsComps(1) ~= metsComps(2) %different compartments, no check for extComp
          if strcmp(metNames{1}, metNames{2}) %the same metabolite, it is a transport reaction
              toIgnoreSimpleTransp(i) = true;
          end
      end
    else %check for advanced 
        metsInd = find(minModel2.S(:,i) ~= 0);
        %check that we have an even number of mets that is larger than 2 and that we don't have a GPR
        if rem(length(metsInd),2) == 0 && length(metsInd) > 2 && isempty(minModel2.grRules{i})
            %Now check that we have pairs of metabolites that are transported
            SVals = full(minModel2.S(metsInd,i));
            metNames = minModel2.metNames(metsInd);
            comps = minModel2.metComps(metsInd);
            success = true;
            while ~isempty(metNames)
                if length(metNames) < 2 %this should not happen
                   error('We should never arrive here!');
                end
                metMatch = find(strcmp(metNames(2:length(metNames)),metNames{1})) + 1;
                if length(metMatch) ~= 1 %we want a match, and only one (if there is some other complex transport rxn, we skip that)
                   success = false;
                   break;
                end
                if (SVals(1) + SVals(metMatch)) ~= 0 %the stoichiometry is not right
                   success = false;
                   break;
                end
                if comps(1) == comps(metMatch) %they must be in different compartments (maybe this need not to be checked)
                   success = false;
                   break;
                end
                %now remove the pair and continue with the next
                SVals([1;metMatch]) = [];
                metNames([1;metMatch]) = [];
                comps([1;metMatch]) = [];
            end
            toIgnoreAdvTransp(i) = success;
        end          
    end
end

%Spontaneous reactions:
toIgnoreSpont = ismember(minModel2.rxns, spontRxnNames);%only 10 rxns in human-GEM

%rxns without GPRs in the s compartment (outside the cell)
sCompInd = find(strcmp(minModel2.comps, extComp));
for i = 1:length(minModel2.rxns)
    metsInd = find(minModel2.S(:,i) ~= 0);
    %check that it has more than 1 metabolite (not exch rxn) and no GPR
    if length(metsInd) > 1 && isempty(minModel2.grRules{i})
        if  all(minModel2.metComps(metsInd) == sCompInd)
            toIgnoreS(i) = true;
        end
    end
end

toIgnoreAllWithoutGPRs = cellfun(@isempty,minModel2.grRules);


%sum(toIgnore)%in total 2375 rxns in human-GEM
%sum(toIgnoreAllTransp)%in total 3337 rxns in human-GEM, so 1,000 more

%Now, try to scale the model to become more favorable for the solver.
%In general, we try to make all fluxes as similar as possible.
%      There is room for improvement here, this function is pretty simple.
%      It only rescales reactions, while it would be possible to rescale
%      metabolites as well (ROS => kROS, albumin => millialbumin, etc., but
%      scaled freely with a mathematical method). It could potentially open up
%      for using less strict margins in the MILP, and make it run faster.
if (skipScaling)
    scaledMinModel = minModel3;
else
    scaledMinModel = rescaleModelForINIT(minModel3);
    scaledMinModel.ub(scaledMinModel.ub > 0) = 1000;
    scaledMinModel.lb(scaledMinModel.lb < 0) = -1000;
end
% create data structure with pre-processing results
prepData.taskReport = taskReport;
prepData.essentialRxns = newEssentialRxns; %essential rxns in the minModel
prepData.taskStruct = taskStruct;
prepData.refModel = cModel;
prepData.minModel = scaledMinModel;
prepData.refModelWithBM = closeModel(cModel); %do this here so we don't have to do this for each sample
prepData.groupIds = groupIds;
prepData.essentialMetsForTasks = essentialMetsForTasks;

prepData.toIgnoreExch = toIgnoreExch;
prepData.toIgnoreImportRxns = toIgnoreImportRxns;
prepData.toIgnoreSimpleTransp = toIgnoreSimpleTransp;
prepData.toIgnoreAdvTransp = toIgnoreAdvTransp;
prepData.toIgnoreSpont = toIgnoreSpont;
prepData.toIgnoreS = toIgnoreS;
prepData.toIgnoreCustomRxns = ismember(minModel2.rxns, customRxnsToIgnore);
prepData.toIgnoreAllWithoutGPRs = toIgnoreAllWithoutGPRs;

end