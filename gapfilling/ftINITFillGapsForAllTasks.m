function [outModel, addedRxns, failedTasks]=ftINITFillGapsForAllTasks(model,refModel,inputFile,printOutput,rxnScores,taskStructure,params,verbose)
% ftINITFillGapsForAllTasks
%   Fills gaps in a model by including reactions from a reference model,
%   so that the resulting model can perform all the tasks in a task list.
%
%   Input:
%   model           model structure
%   refModel        reference model from which to include reactions. It is
%                   expected to already contain the reactions in model, so
%                   that no merge is needed for each task
%   inputFile       a task list in Excel format. See the function
%                   parseTaskList for details (optional if taskStructure is
%                   supplied)
%   printOutput     true if the results of the test should be displayed
%                   (optional, default true)
%   rxnScores       scores for each of the reactions in the reference
%                   model. Only negative scores are allowed. The solver will
%                   try to maximize the sum of the scores for the included
%                   reactions (optional, default is -1 for all reactions)
%   taskStructure   structure with the tasks, as from parseTaskList. If
%                   this is supplied then inputFile is ignored
%   params          *obsolete option*
%   verbose         if true, the MILP progression will be shown (optional,
%                   default false)
%
%
%   Output:
%   outModel        model structure with reactions added to perform the
%                   tasks
%   addedRxns       MxN matrix with the added reactions (M) from refModel
%                   for each task (N). An element is true if the corresponding
%                   reaction is added in the corresponding task.
%                   Failed tasks and SHOULD FAIL tasks are ignored
%   failedTasks     Nx1 logical, true for each task that could not be
%                   gap-filled, either because no feasible solution exists
%                   using the reference model or because the attempt threw.
%                   Such tasks add no reactions, so they are otherwise
%                   indistinguishable from tasks that already worked
%
%   This is fitTasks with gapFillMode set to 'preMerged', which is the mode
%   that suits ftINIT: the task constraints are written into the reference
%   model as well, and ftINITFillGaps replaces fillGaps.
%
%   See also: fitTasks, which this function wraps.
%
% Usage: [outModel, addedRxns, failedTasks]=ftINITFillGapsForAllTasks(model,...
%           refModel,inputFile,printOutput,rxnScores,taskStructure,params,verbose)

if nargin<4 || isempty(printOutput)
    printOutput=true;
end
if nargin<5
    rxnScores=[];
end
if nargin<6
    taskStructure=[];
end
if nargin<7
    params=[];
end
if nargin<8 || isempty(verbose)
    verbose=false;
end

[outModel, addedRxns, failedTasks]=fitTasks(model,refModel,inputFile, ...
    'printOutput',printOutput, ...
    'rxnScores',rxnScores, ...
    'taskStructure',taskStructure, ...
    'gapFillMode','preMerged', ...
    'params',params, ...
    'verbose',verbose);
end
