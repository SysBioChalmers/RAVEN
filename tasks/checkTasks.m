function [taskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,inputFile,varargin)
% checkTasks  Perform a set of simulations as defined in a task file.
%
% This function is used for defining a set of tasks for a model to perform.
% The tasks are defined by defining constraints on the model, and if the
% problem is feasible, then the task is considered successful. In general,
% each row can contain one constraint on uptakes, one constraint on outputs,
% one new equation, and one change of reaction bounds. If more bounds are
% needed to define the task, then several rows can be used for each task.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% inputFile : char
%     a task list in Excel format. See the function parseTaskList for
%     details (optional if taskStructure is supplied).
%
% Name-Value Arguments
% --------------------
% printOutput : logical
%     true if the results of the test should be displayed (default true).
% printOnlyFailed : logical
%     true if only tasks that failed should be displayed (default false).
% getEssential : logical
%     true if the essential reactions should be calculated for all the
%     tasks. This option is used with runINIT (default false).
% taskStructure : struct
%     structure with the tasks, as from parseTaskList. If this is supplied
%     then inputFile is ignored.
% runParallel : logical
%     true to evaluate tasks in parallel workers, since each task is
%     independent of the others (see parallelWorkersRAVEN). Gurobi is
%     automatically pinned to one thread per worker by optimizeProb, so
%     this does not oversubscribe cores. Default false, matching the
%     previous (serial) behaviour, since checkTasks is often called on
%     small task lists where starting a pool would only add overhead.
%
% Returns
% -------
% taskReport : struct
%     structure with the results, with fields:
%
%     - id : cell array with the id of the task
%     - description : cell array with the description of the task
%     - ok : boolean array with true if the task was successful
% essentialRxns : logical
%     MxN matrix with the essential reactions (M) for each task (N). An
%     element is true if the corresponding reaction is essential in the
%     corresponding task. Failed tasks and SHOULD FAIL tasks are ignored.
%     This is used by the INIT algorithm (if tasks are supplied). If
%     getEssential=false then essentialRxns=false(nRxns,nTasks).
% taskStructure : struct
%     structure with the tasks, as from parseTaskList.
% essentialFluxes : double
%     the fluxes of the essential rxns - same structure as essentialRxns.
%
% Examples
% --------
%     [taskReport, essentialRxns, taskStructure] = checkTasks(model, inputFile, ...
%         printOutput, printOnlyFailed, getEssential, taskStructure);

p=parseRAVENargs(varargin, {'printOutput',true; 'printOnlyFailed',false; 'getEssential',false; 'taskStructure',[]; 'runParallel',false});
printOutput=p.printOutput;
printOnlyFailed=p.printOnlyFailed;
getEssential=p.getEssential;
taskStructure=p.taskStructure;
runParallel=p.runParallel;

%Prepare the input model a little
model.b=zeros(numel(model.mets),2);

modelMets=upper(strcat(model.metNames,'[',model.comps(model.metComps),']'));
if ~isfield(model,'unconstrained')
    EM='Exchange metabolites should normally not be removed from the model when using checkTasks. Inputs and outputs are defined in the task file instead. Use importModel(file,false) to import a model with exchange metabolites remaining';
    warning('RAVEN:warning', '%s', EM);
end

%Parse the task file
if isempty(taskStructure)
    taskStructure=parseTaskList(inputFile);
end

nTasks=numel(taskStructure);
essentialRxns=false(numel(model.rxns),nTasks);
essentialFluxes = NaN(numel(model.rxns),nTasks);

%Each task is independent (tModel is rebuilt from model at the top of
%every iteration, nothing carries over), so this loop can run in
%parallel. ids/descriptions/ok are plain arrays sliced by i and only
%assembled into taskReport once the loop is done, because parfor's
%classifier rejects slicing a struct field that is itself indexed with
%cell/array subscripts (taskReport.id{i,1}=... is not a supported sliced
%form). runParallel=false (the default) makes this run exactly like the
%previous for loop, in the client, serially.
ids=cell(nTasks,1);
descriptions=cell(nTasks,1);
ok=false(nTasks,1);

nW=parallelWorkersRAVEN(runParallel);
parfor (i=1:nTasks, nW)
    tModel=model;
    ids{i}=taskStructure(i).id;
    descriptions{i}=taskStructure(i).description;
    %Set the inputs
    if ~isempty(taskStructure(i).inputs)
        [I, J]=ismember(upper(taskStructure(i).inputs),modelMets);
        J=J(I); %Only keep the ones with matches
        K=ismember(upper(taskStructure(i).inputs),'ALLMETS');
        L=~cellfun('isempty',strfind(upper(taskStructure(i).inputs),'ALLMETSIN'));
        %Check that all metabolites are either real metabolites or
        %ALLMETS/ALLMETSIN
        if ~all(I|K|L)
            fprintf(['ERROR: Could not find all inputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
            ok(i)=false;
            continue;
        end
        if numel(J)~=numel(unique(J))
            EM=['The constraints on some input(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
            error('RAVEN:badInput', '%s', EM);
        end
        %If all metabolites should be added
        if any(K)
            %Check if ALLMETS is the first metabolite. Otherwise print a
            %warning since it will write over any other constraints that
            %are set
            if K(1)==0
                EM=['ALLMETS is used as an input in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                warning('RAVEN:warning', '%s', EM);
            end
            %Use the first match of ALLMETS. There should only be one, but
            %still..
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
                    error('RAVEN:badInput', '%s', EM);
                end
            end
        end
        %Then add the normal constraints
        if any(J)
            tModel.b(J,1)=taskStructure(i).UBin(I)*-1;
            tModel.b(J,2)=taskStructure(i).LBin(I)*-1;
        end
    end
    %Set the outputs
    if ~isempty(taskStructure(i).outputs)
        [I, J]=ismember(upper(taskStructure(i).outputs),modelMets);
        J=J(I); %Only keep the ones with matches
        K=ismember(upper(taskStructure(i).outputs),'ALLMETS');
        L=~cellfun('isempty',strfind(upper(taskStructure(i).outputs),'ALLMETSIN'));
        %Check that all metabolites are either real metabolites or
        %ALLMETS/ALLMETSIN
        if ~all(I|K|L)
            fprintf(['ERROR: Could not find all outputs in "[' taskStructure(i).id '] ' taskStructure(i).description '"\n']);
            ok(i)=false;
            continue;
        end
        if numel(J)~=numel(unique(J))
            EM=['The constraints on some output(s) in "[' taskStructure(i).id '] ' taskStructure(i).description '" are defined more than one time'];
            error('RAVEN:badInput', '%s', EM);
        end
        %If all metabolites should be added
        if any(K)
            %Check if ALLMETS is the first metabolite. Otherwise print a
            %warning since it will write over any other constraints that
            %are set
            if K(1)==0
                EM=['ALLMETS is used as an output in "[' taskStructure(i).id '] ' taskStructure(i).description '" but it it not the first metabolite in the list. Constraints defined for the metabolites before it will be over-written'];
                warning('RAVEN:warning', '%s', EM);
            end
            %Use the first match of ALLMETS. There should only be one, but
            %still..
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
                    error('RAVEN:badInput', '%s', EM);
                end
            end
        end
        %Then add the normal constraints
        if any(J)
            %Verify that IN and OUT bounds are consistent. Cannot require
            %that a metabolite is simultaneously input AND output at some
            %nonzero flux.
            I = find(I);  % otherwise indexing becomes confusing
            nonzero_LBin = tModel.b(J,2) < 0;
            nonzero_LBout = taskStructure(i).LBout(I) > 0;
            if any(nonzero_LBin & nonzero_LBout)
                EM=['The IN LB and OUT LB in "[' taskStructure(i).id '] ' taskStructure(i).description '" cannot be nonzero for the same metabolite'];
                error('RAVEN:badInput', '%s', EM);
            end
            tModel.b(J(nonzero_LBout),1)=taskStructure(i).LBout(I(nonzero_LBout));
            tModel.b(J,2)=taskStructure(i).UBout(I);
        end
    end
    %Add new rxns
    if ~isempty(taskStructure(i).equations)
        % num2str on the whole column right-aligns every row to a common
        % width, embedding a leading space in "TEMPORARY_ 1" once any id
        % reaches two digits ("TEMPORARY_10"); format each one independently.
        rxnIds=arrayfun(@(x) sprintf('TEMPORARY_%d',x), ...
            (1:numel(taskStructure(i).equations))', 'UniformOutput', false);
        % Built as a single struct() call, not incremental rxn.field=...
        % assignments: parfor's variable classifier cannot establish that a
        % variable built up via several dot-indexed assignments is a
        % same-shape temporary in every iteration, and rejects it outright.
        rxn=struct('equations',{taskStructure(i).equations}, ...
            'lb',taskStructure(i).LBequ, 'ub',taskStructure(i).UBequ, ...
            'rxns',{rxnIds});
        %Allow for new metabolites to be added. This is because it should
        %be possible to add, say, a whole new pathway
        tModel=addRxns(tModel,rxn,3,[],true);
    end
    %Add changed bounds
    if ~isempty(taskStructure(i).changed)
        tModel=setParam(tModel,'lb',taskStructure(i).changed,taskStructure(i).LBrxn);
        tModel=setParam(tModel,'ub',taskStructure(i).changed,taskStructure(i).UBrxn);
    end
    
    %Solve and print
    sol=solveLP(tModel);
    if ~isempty(sol.x)
        %assign the fluxes
        essentialFluxes(:,i) = sol.x(1:numel(model.rxns));
        
        if ~taskStructure(i).shouldFail
            ok(i)=true;
            if printOnlyFailed==false && printOutput==true
                fprintf(['PASS: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
            %Calculate the essential reactions
            if getEssential==true
                [~, taskEssential]=getEssentialRxns(tModel);
                %This is because there could be more reactions in tModel
                %than in model. Built as a full local column and written
                %with a plain ':' row index: parfor only accepts a sliced
                %output write when every subscript besides the loop
                %variable is ':', not a data-dependent index vector.
                essentialCol=false(numel(model.rxns),1);
                essentialCol(taskEssential(taskEssential<=numel(model.rxns)))=true;
                essentialRxns(:,i)=essentialCol;
            end
        else
            ok(i)=false;
            if printOutput==true
                fprintf(['PASS (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        end
    else
        if ~taskStructure(i).shouldFail
            ok(i)=false;
            if printOutput==true
                fprintf(['FAIL: [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        else
            ok(i)=true;
            if printOnlyFailed==false && printOutput==true
                fprintf(['FAIL (should fail): [' taskStructure(i).id '] ' taskStructure(i).description '\n']);
            end
        end
    end
    if taskStructure(i).printFluxes && ~isempty(sol.x)
        sol=solveLP(tModel,1);
        if ~isempty(sol.x)
            printFluxes(tModel,sol.x,false,10^-6,[],'%rxnID (%eqn):%flux\n');
            fprintf('\n');
        end
    end
end

taskReport.id=ids;
taskReport.description=descriptions;
taskReport.ok=ok;

end
