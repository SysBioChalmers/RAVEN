function [outModel, addedRxns]=fitTasks(model,refModel,inputFile,printOutput,rxnScores,taskStructure,params)
% fitTasks
%   Fills gaps in a model by including reactions from a reference model,
%   so that the resulting model can perform all the tasks in a task list
%
%   Input:
%   model           model structure
%   refModel        reference model from which to include reactions
%   inputFile       a task list in Excel format. See the function
%                   parseTaskList for details (opt if taskStructure is
%                   supplied)
%   printOutput     true if the results of the test should be displayed
%                   (opt, default true)
%   rxnScores       scores for each of the reactions in the reference
%                   model. Only negative scores are allowed. The solver will
%                   try to maximize the sum of the scores for the included
%                   reactions (opt, default is -1 for all reactions)
%   taskStructure   structure with the tasks, as from parseTaskList. If
%                   this is supplied then inputFile is ignored (opt)
%   params          parameter structure as used by getMILPParams (opt)
%
%
%   Output:
%   outModel        model structure with reactions added to perform the
%                   tasks
%   addedRxns       MxN matrix with the added reactions (M) from refModel
%                   for each task (N). An element is true if the corresponding
%                   reaction is added in the corresponding task.
%                   Failed tasks and SHOULD FAIL tasks are ignored
%
%   This function fills gaps in a model by using a reference model, so
%   that the resulting model can perform a list of metabolic tasks. The
%   gap-filling is done in a task-by-task manner, rather than solving for
%   all tasks at once. This means that the order of the tasks could influence
%   the result.
%
%   Usage: [outModel, addedRxns]=fitTasks(model,refModel,inputFile,printOutput,...
%           rxnScores,taskStructure,params)

if nargin<4
    printOutput=true;
end
if nargin<5
    rxnScores=ones(numel(refModel.rxns),1)*-1;
end
if isempty(rxnScores)
    rxnScores=ones(numel(refModel.rxns),1)*-1;
end
if nargin<6
    taskStructure=[];
end
if nargin<7
    params=[];
end

if isempty(taskStructure) && ~isfile(inputFile)
    error('Task file %s cannot be found',string(inputFile));
end

if strcmpi(model.id,refModel.id)
    fprintf('NOTE: The model and reference model have the same IDs. The ID for the reference model was set to "refModel" in order to keep track of the origin of reactions.\n');
    refModel.id='refModel';
end

if any(rxnScores>=0)
    EM='Only negative values are allowed in rxnScores';
    dispEM(EM);
end

%Prepare the input models a little
model.b=zeros(numel(model.mets),2);
modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
%This is the mets in the reference model. Used if the tasks involve
%metabolites that doesn't exist in the model
largeModelMets=upper(strcat(refModel.metNames,'[',refModel.comps(refModel.metComps),']'));

if ~isfield(model,'unconstrained')
    EM='Exchange metabolites should normally not be removed from the model when using checkTasks. Inputs and outputs are defined in the task file instead. Use importModel(file,false) to import a model with exchange metabolites remaining';
    dispEM(EM,false);
end

if isempty(taskStructure)
    taskStructure=parseTaskList(inputFile);
end

tModel=model;
addedRxns=false(numel(refModel.rxns),numel(taskStructure));
supressWarnings=false;
nAdded=0;
for i=1:numel(taskStructure)
    if ~taskStructure(i).shouldFail
        %Set the inputs
        if ~isempty(taskStructure(i).inputs)
            [I, J]=ismember(upper(taskStructure(i).inputs),modelMets);
            K=ismember(upper(taskStructure(i).inputs),'ALLMETS');
            L=~cellfun('isempty',strfind(upper(taskStructure(i).inputs),'ALLMETSIN'));
            %Check that all metabolites are either real metabolites or
            %ALLMETS/ALLMETSIN
            goodMets=I|K|L;
            if ~all(goodMets)
                %Not all of the inputs could be found in the small model.
                %Check if they exist in the large model
                [found, metMatch]=ismember(upper(taskStructure(i).inputs(~goodMets)),largeModelMets);
                if ~all(found)
                    EM=['Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '" in either model'];
                    disp(EM);
                else
                    %Otherwise add them to the model
                    met.metNames=refModel.metNames(metMatch);
                    met.compartments=refModel.comps(refModel.metComps(metMatch));
                    
                    %Add the metabolite both to the base model and the
                    %model used in the current task
                    model=addMets(model,met);
                    tModel=addMets(tModel,met);
                    modelMets=[modelMets;upper(taskStructure(i).inputs(~goodMets))];
                end
                
                %By now the indexes might be getting a bit confusing, but
                %this is to update the indexes of the "real" metabolites to
                %point to the newly added ones
                I(~goodMets)=true; %All the bad ones are fixed at this stage
                J(~goodMets)=numel(modelMets)-numel(metMatch)+1:numel(modelMets);
            end
            if numel(J(I))~=numel(unique(J(I)))
                EM=['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
                dispEM(EM);
            end
            %If all metabolites should be added
            if any(K)
                %Check if ALLMETS is the first metabolite. Otherwise print
                %a warning since it will write over any other constraints
                %that are set
                if K(1)==0
                    EM=['ALLMETS is used as an input in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                    dispEM(EM,false);
                end
                %Use the first match of ALLMETS. There should only be one,
                %but still..
                tModel.b(:,1)=taskStructure(i).UBin(find(K,1))*-1;
            end
            %If metabolites in a specific compartment should be used
            if any(L)
                L=find(L);
                for j=1:numel(L)
                    %The compartment defined
                    compartment=upper(taskStructure(i).inputs{L(j)}(11:end-1));
                    %Check if it exists in the model
                    C=find(ismember(upper(model.comps),compartment));
                    if any(C)
                        %Match to metabolites
                        tModel.b(model.metComps==C,1)=taskStructure(i).UBin(L(j))*-1;
                    else
                        EM=['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist'];
                        dispEM(EM);
                    end
                end
            end
            %Then add the normal constraints
            if any(J(I))
                tModel.b(J(I),1)=taskStructure(i).UBin(I)*-1;
                tModel.b(J(I),2)=taskStructure(i).LBin(I)*-1;
            end
        end
        %Set the outputs
        if ~isempty(taskStructure(i).outputs)
            [I, J]=ismember(upper(taskStructure(i).outputs),modelMets);
            K=ismember(upper(taskStructure(i).outputs),'ALLMETS');
            L=~cellfun('isempty',strfind(upper(taskStructure(i).outputs),'ALLMETSIN'));
            %Check that all metabolites are either real metabolites or
            %ALLMETS/ALLMETSIN
            goodMets=I|K|L;
            if ~all(goodMets)
                %Not all of the outputs could be found in the small model.
                %Check if they exist in the large model
                [found, metMatch]=ismember(upper(taskStructure(i).outputs(~goodMets)),largeModelMets);
                if ~all(found)
                    EM=['Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '" in either model'];
                    dispEM(EM);
                else
                    %Otherwise add them to the model
                    met.metNames=refModel.metNames(metMatch);
                    met.compartments=refModel.comps(refModel.metComps(metMatch));
                    
                    %Add the metabolite both to the base model and the
                    %model used in the current task
                    model=addMets(model,met);
                    tModel=addMets(tModel,met);
                    modelMets=[modelMets;upper(taskStructure(i).outputs(~goodMets))];
                end
                
                %By now the indexes might be getting a bit confusing, but
                %this is to update the indexes of the "real" metabolites to
                %point to the newly added ones
                I(~goodMets)=true; %All the bad ones are fixed at this stage
                J(~goodMets)=numel(modelMets)-numel(metMatch)+1:numel(modelMets);
            end
            if numel(J(I))~=numel(unique(J(I)))
                EM=['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
                dispEM(EM);
            end
            %If all metabolites should be added
            if any(K)
                %Check if ALLMETS is the first metabolite. Otherwise print
                %a warning since it will write over any other constraints
                %that are set
                if K(1)==0
                    EM=['ALLMETS is used as an output in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                    dispEM(EM,false);
                end
                %Use the first match of ALLMETS. There should only be one,
                %but still..
                tModel.b(:,2)=taskStructure(i).UBout(find(K,1));
            end
            %If metabolites in a specific compartment should be used
            if any(L)
                L=find(L);
                for j=1:numel(L)
                    %The compartment defined
                    compartment=upper(taskStructure(i).outputs{L(j)}(11:end-1));
                    %Check if it exists in the model
                    C=find(ismember(upper(model.comps),compartment));
                    if any(C)
                        %Match to metabolites
                        tModel.b(model.metComps==C,2)=taskStructure(i).UBout(L(j));
                    else
                        EM=['The compartment defined for ALLMETSIN in "[' taskStructure(i).id '] ' taskStructure(i).description '" does not exist'];
                        dispEM(EM);
                    end
                end
            end
            %Then add the normal constraints
            if any(J(I))
                %Verify that IN and OUT bounds are consistent. Cannot require
                %that a metabolite is simultaneously input AND output at some
                %nonzero flux.
                J = J(I);
                I = find(I);  % otherwise indexing becomes confusing
                nonzero_LBin = tModel.b(J,2) < 0;
                nonzero_LBout = taskStructure(i).LBout(I) > 0;
                if any(nonzero_LBin & nonzero_LBout)
                    EM=['The IN LB and OUT LB in "[' taskStructure(i).id '] ' taskStructure(i).description '" cannot be nonzero for the same metabolite'];
                    dispEM(EM);
                end
                tModel.b(J(nonzero_LBout),1)=taskStructure(i).LBout(I(nonzero_LBout));
                tModel.b(J,2)=taskStructure(i).UBout(I);
            end
        end
        
        %Add new rxns
        if ~isempty(taskStructure(i).equations)
            rxn.equations=taskStructure(i).equations;
            rxn.lb=taskStructure(i).LBequ;
            rxn.ub=taskStructure(i).UBequ;
            rxn.rxns=strcat({'TEMPORARY_'},num2str((1:numel(taskStructure(i).equations))'));
            tModel=addRxns(tModel,rxn,3);
        end
        %Add changed bounds
        if ~isempty(taskStructure(i).changed)
            tModel=setParam(tModel,'lb',taskStructure(i).changed,taskStructure(i).LBrxn);
            tModel=setParam(tModel,'ub',taskStructure(i).changed,taskStructure(i).UBrxn);
        end
        
        %Solve and print. Display a warning if the problem is not solveable
        sol=solveLP(tModel);
        if isempty(sol.x)
            %Only do gap-filling if it cannot be solved
            failed=false;
            try
                [~, ~, newRxns, newModel, exitFlag]=fillGaps(tModel,refModel,false,true,supressWarnings,rxnScores,params);
                if exitFlag==-2
                    EM=['"[' taskStructure(i).id '] ' taskStructure(i).description '" was aborted before reaching optimality. Consider increasing params.maxTime\n'];
                    dispEM(EM,false);
                end
            catch
                EM=['"[' taskStructure(i).id '] ' taskStructure(i).description '" could not be performed for any set of reactions\n'];
                dispEM(EM,false);
                failed=true;
            end
            if failed==false
                if ~isempty(newRxns)
                    nAdded=nAdded+numel(newRxns);
                    
                    %Add the reactions to the base model. It is not correct
                    %to use newModel directly, as it may contain
                    %reactions/constraints that are specific to this task
                    model=mergeModels({model,removeReactions(newModel,setdiff(newModel.rxns,newRxns),true,true)},'metNames',true);
                    
                    %Keep track of the added reactions
                    addedRxns(ismember(refModel.rxns,newRxns),i)=true;
                end
                if printOutput==true
                    fprintf(['[' taskStructure(i).id '] ' taskStructure(i).description ': Added ' num2str(numel(newRxns)) ' reaction(s), ' num2str(nAdded) ' reactions added in total\n']);
                end
            end
        else
            if printOutput==true
                fprintf(['[' taskStructure(i).id '] ' taskStructure(i).description ': Added 0 reaction(s), ' num2str(nAdded) ' reactions added in total\n']);
            end
        end
        supressWarnings=true;
        
        %Print the output if chosen
        if taskStructure(i).printFluxes && printOutput
            if ~isempty(sol.x)
                sol=solveLP(tModel,1);
                printFluxes(tModel,sol.x,false,10^-5,[],'%rxnID (%eqn):%flux\n');
                fprintf('\n');
            else
                %If the problem wasn't solveable then the gap-filled model
                %should be used
                if failed==false
                    sol=solveLP(newModel,1);
                    printFluxes(newModel,sol.x,false,10^-5,[],'%rxnID (%eqn):%flux\n');
                    fprintf('\n');
                end
            end
        end
        
        tModel=model;
        %Since new mets are added by merging the new reactions and not only
        %from the task sheet
        modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
    else
        EM=['"[' taskStructure(i).id '] ' taskStructure(i).description '" is set as SHOULD FAIL. Such tasks cannot be modelled using this approach and the task is therefore ignored\n'];
        dispEM(EM,false);
    end
end
model.b(:,2) = [];  % resume field b
outModel=model;
end
