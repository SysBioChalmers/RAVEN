%run this test case with the command
%results = runtests('checkTasksTests.m')
function tests = checkTasksTests
    tests = functiontests(localfunctions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModel = getTaskTstModel()
    testModel = struct();
    testModel.id = 'testModel';
    testModel.rxns = {};
    testModel.S=[];
    testModel.rev=[];
    testModel.mets = {'as';'ac';'bc';'cc';'dc';'ec';'es';'fc'};
    testModel.metNames = {'a';'a';'b';'c';'d';'e';'e';'f'};
    testModel.comps = {'s';'c'};
    testModel.compNames = testModel.comps;
    testModel.metComps = [1;2;2;2;2;2;1;2];
    testModel.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
    testModel.grRules = {};
    testModel.rxnGeneMat = [];

    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
    rxnsToAdd.equations = {'=> a[s]';...
                           'a[s] <=> a[c]';...
                           'a[c] <=> b[c] + c[c]';...
                           'a[c] <=> 2 d[c]';...
                           'b[c] + c[c] => e[c]';...
                           '2 d[c] => e[c]';...
                           'e[c] => e[s]';...
                           'e[s] =>';...
                           'a[c] <=> f[c]';...
                           'f[c] <=> e[c]'};
    rxnsToAdd.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};
    testModel = addRxns(testModel,rxnsToAdd, 3);
    testModel.c = [0;0;0;0;0;0;0;1;0;0];%optimize for output flux, if this is used, not sure
    testModel.ub = repmat(1000,10,1);
    testModel.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
    testModel.rxnNames = testModel.rxns;
    testModel.b = repmat(0,8,1);
end

function testModelTasks = getTaskTstModelTasks()
    testModelTasks = struct();
    testModelTasks.id = 'Gen e[s] from a[s]';
    testModelTasks.description = 'Gen e[s] from a[s]';
    testModelTasks.shouldFail = false;
    testModelTasks.printFluxes = false;
    testModelTasks.comments = '';
    testModelTasks.inputs = {'a[s]'};
    testModelTasks.LBin = 0;
    testModelTasks.UBin = inf;
    testModelTasks.outputs = {'e[s]'};
    testModelTasks.LBout = 1;
    testModelTasks.UBout = 1;
    testModelTasks.equations = {};
    testModelTasks.LBequ = [];
    testModelTasks.UBequ = [];
    testModelTasks.changed = {};
    testModelTasks.LBrxn = {};
    testModelTasks.UBrxn = {};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TXXXX: checkTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testcheckTasks(testCase)

    testModel = getTaskTstModel();
    testModel = closeModel(testModel);
    testModelTasks = getTaskTstModelTasks();
    
    %the task should work with the full model
    taskReport = checkTasks(testModel,[],true,false,false,testModelTasks);
    verifyTrue(testCase, taskReport.ok == 1)
    
    %now remove R2, then it should fail
    testModel2 = removeReactions(testModel, {'R2'});
    taskReport2 = checkTasks(testModel2,[],true,false,false,testModelTasks);
    verifyTrue(testCase, taskReport2.ok == 0)
end


