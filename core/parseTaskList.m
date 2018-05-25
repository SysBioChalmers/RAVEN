function taskStruct=parseTaskList(inputFile)
% parseTaskList
%   Parses a task list file.
%
%   inputFile       a task list in Excel format. The file must contain a
%                   sheet named TASKS, which in turn may contain the
%                   following column headers (note, all rows starting with
%                   a non-empty cell are removed. The first row after that
%                   is considered the headers):
%                   ID
%                       the only required header. Each task must have a
%                       unique id (string or numeric). Tasks can span multiple
%                       rows, only the first row in each task should have
%                       an id
%                   DESCRIPTION
%                       description of the task
%                   IN
%                       allowed input(s) for the task. Metabolite names
%                       should be on the form
%                       "model.metName[model.comps]". Several inputs
%                       can be delimited by ";". If so, then the same
%                       bounds are used for all inputs. If that is not
%                       wanted, then use several rows for the task
%                   IN LB
%                       lower bound for the uptake of the metabolites in
%                       the row (opt, default 0 which corresponds to a
%                       minimal uptake of 0 units)
%                   IN UB
%                       upper bound for the uptake of the metabolites in
%                       the row (opt, default 1000 which corresponds to a
%                       maximal uptake of 1000 units)
%                   OUT
%                       allowed output(s) for the task (see IN)
%                   OUT LB
%                       lower bound for the production of the metabolites in
%                       the row (opt, default 0 which corresponds to a
%                       minimal production of 0 units)
%                   OUT UB
%                       upper bound for the production of the metabolites in
%                       the row (opt, default 1000 which corresponds to a
%                       maximal production of 1000 units)
%                   EQU
%                       equation to add. The equation should be on the form
%                       "0.4 A + 2 B <=> (or =>) C" and the metabolites
%                       should be on the form
%                       "model.metName[model.comps]" (opt)
%                   EQU LB
%                       lower bound for the equation (opt, default -1000
%                       for reversible and 0 for irreversible)
%                   EQU UB
%                       upper bound for the equation (opt, default 1000)
%                   CHANGED RXN
%                       reaction ID for which to change the bounds for.
%                       Several IDs can be delimited by ";". If so,
%                       then the same bounds are used for all reactions. If
%                       that is not wanted, then use several rows for the task
%                   CHANGED LB
%                       lower bound for the reaction
%                   CHANGED UB
%                       upper bound for the reaction
%                   SHOULD FAIL
%                       true if the correct behavior of the model is to
%                       not have a feasible solution given the constraints
%                       (opt, default false)
%                   PRINT FLUX
%                       true if the function should print the corresponding
%                       flux distribution for a task. Can be useful for
%                       testing (opt, default false)
%
%   taskStruct      array of structures with the following fields
%       id          the id of the task
%       description the description of the task
%       shouldFail  true if the task should fail
%       printFluxes true if the fluxes should be printed
%       comments    string with comments
%       inputs      cell array with input metabolites (in the form metName[comps])
%       LBin        array with lower bounds on inputs (default, 0)
%       UBin        array with upper bounds on inputs (default, 1000)
%       outputs     cell array with output metabolites (in the form metName[comps])
%       LBout       array with lower bounds on outputs (default, 0)
%       UBout       array with upper bounds on outputs (default, 1000)
%       equations   cell array with equations (with mets in the form metName[comps])
%       LBequ       array with lower bounds on equations (default, -1000 for
%                   reversible and 0 for irreversible)
%       UBequ       array with upper bounds on equations (default, 1000)
%       changed     cell array with reactions to change bounds for
%       LBrxn       array with lower bounds on changed reactions
%       UBrxn       array with upper bounds on changed reactions
%
%   This function is used for defining a set of tasks for a model to
%   perform. The tasks are defined by defining constraints on the model,
%   and if the problem is feasible, then the task is considered successful.
%   In general, each row can contain one constraint on uptakes, one
%   constraint on outputs, one new equation, and one change of reaction
%   bounds. If more bounds are needed to define the task, then several rows
%   can be used for each task. To perform the task use checkTasks or
%   fitTasks.
%
%   NOTE: The general metabolites "ALLMETS" and "ALLMETSIN[comps]"
%   can be used as inputs or outputs in the similar manner to normal
%   metabolites. This is a convenient way to, for example, allow excretion of
%   all metabolites to check whether it's the synthesis of some metabolite
%   that is limiting or whether it's the degradation of some byproduct. One
%   important difference is that only the upper bounds are used for these general
%   metabolites. That is, you can only say that uptake or excretion is
%   allowed, not that it is required. This is to avoid conflicts where the
%   constraints for the general metabolites overwrite those of the real
%   ones.
%
%   Usage: taskStruct=parseTaskList(inputFile)
%
%   Rasmus Agren, 2017-02-28
%

%Load the tasks file
[raw,flag]=loadSheet(loadWorkbook(inputFile), 'TASKS');
if flag~=0
    EM=['Could not load sheet "TASKS" from ' inputFile];
    dispEM(EM);
end

%Remove all lines starting with "#" (or actually any character) and all
%empty columns
raw=cleanSheet(raw);

%Captions
columns={'ID';'DESCRIPTION';'IN';'IN LB';'IN UB';'OUT';'OUT LB';'OUT UB';'EQU';'EQU LB';'EQU UB';'CHANGED RXN';'CHANGED LB';'CHANGED UB';'SHOULD FAIL';'PRINT FLUX';'COMMENTS'};

%Match the columns
[I, colI]=ismember(columns,raw(1,:));

%Check that the ID field is present
if I(1)==0
    EM='The TASKS sheet must have a column named ID';
    dispEM(EM);
end

%Add default bounds where needed
for i=[4 5 7 8]
    I=cellfun(@isempty,raw(:,colI(i)));
    if i==5 || i==8
        raw(I,colI(i))={1000};
    else
        raw(I,colI(i))={0};
    end
end

%Create an empty task structure
eTask.id='';
eTask.description='';
eTask.shouldFail=false;
eTask.printFluxes=false;
eTask.comments='';
eTask.inputs={};
eTask.LBin=[];
eTask.UBin=[];
eTask.outputs={};
eTask.LBout=[];
eTask.UBout=[];
eTask.equations={};
eTask.LBequ=[];
eTask.UBequ=[];
eTask.changed={};
eTask.LBrxn=[];
eTask.UBrxn=[];

%Main loop
taskStruct=[];
task=eTask;
if isnumeric(raw{2,colI(1)})
    task.id=num2str(raw{2,colI(1)});
else
    task.id=raw{2,colI(1)};
end
task.description=raw{2,colI(2)};
if ~isempty(raw{2,colI(15)})
    task.shouldFail=true;
end
if ~isempty(raw{2,colI(16)})
    task.printFluxes=true;
end
if ~isempty(raw{2,colI(17)})
    task.comments=raw{2,colI(17)};
end

for i=2:size(raw,1)
    %Set the inputs
    if ischar(raw{i,colI(3)})
        inputs=regexp(raw{i,colI(3)},';','split');
        task.inputs=[task.inputs;inputs(:)];
        task.LBin=[task.LBin;ones(numel(inputs),1)*raw{i,colI(4)}];
        task.UBin=[task.UBin;ones(numel(inputs),1)*raw{i,colI(5)}];
    end
    %Set the outputs
    if ischar(raw{i,colI(6)})
        outputs=regexp(raw{i,colI(6)},';','split');
        task.outputs=[task.outputs;outputs(:)];
        task.LBout=[task.LBout;ones(numel(outputs),1)*raw{i,colI(7)}];
        task.UBout=[task.UBout;ones(numel(outputs),1)*raw{i,colI(8)}];
    end
    %Add new rxns
    if ischar(raw{i,colI(9)})
        task.equations=[task.equations;raw{i,colI(9)}];
        if ~isempty(raw{i,colI(10)})
            task.LBequ=[task.LBequ;raw{i,colI(10)}];
        else
            if any(strfind(raw{i,colI(9)},'<=>'))
                task.LBequ=[task.LBequ;-1000];
            else
                task.LBequ=[task.LBequ;0];
            end
        end
        if ~isempty(raw{i,colI(11)})
            task.UBequ=[task.UBequ;raw{i,colI(11)}];
        else
            task.UBequ=[task.UBequ;1000];
        end
    end
    %Add changed bounds
    if ischar(raw{i,colI(12)})
    	changed=regexp(raw{i,colI(12)},';','split');
        task.changed=[task.changed;changed(:)];
        task.LBrxn=[task.LBrxn;ones(numel(changed),1)*raw{i,colI(13)}];
        task.UBrxn=[task.UBrxn;ones(numel(changed),1)*raw{i,colI(14)}];
    end

    %Check if it should add more constraints
    if i<size(raw,1)
        if isempty(raw{i+1,colI(1)})
            continue;
        end
    end

    taskStruct=[taskStruct;task];
    task=eTask;
    if i<size(raw,1)
        if isnumeric(raw{i+1,colI(1)})
            task.id=num2str(raw{i+1,colI(1)});
        else
            task.id=raw{i+1,colI(1)};
        end
        task.description=raw{i+1,colI(2)};
        if ~isempty(raw{i+1,colI(15)})
            task.shouldFail=true;
        end
        if ~isempty(raw{i+1,colI(16)})
            task.printFluxes=true;
        end
        if ~isempty(raw{i+1,colI(17)})
            task.comments=raw{i+1,colI(17)};
        end
    end
end

%Should add more checks, such as unique IDs and missing headers

end
