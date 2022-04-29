%run this test case with the command
%results = runtests('tinitTests.m')
function tests = tinitTests
    tests = functiontests(localfunctions);
end

function testparsexpTask1List(testCase)
    sourceDir = fileparts(which(mfilename));
    taskStruct = parseTaskList(strcat(sourceDir, '/test_data/test_tasks.txt'));
    taskStructExcel = parseTaskList(strcat(sourceDir, '/test_data/test_tasks.xls'));
    %check that all fields in the first line are what we expect
    
    expTask1.id='ER';
    expTask1.description='Aerobic rephosphorylation of ATP from glucose';
    expTask1.shouldFail=true;
    expTask1.printFluxes=true;
    expTask1.comments='Messed up reaction';
    expTask1.inputs={'O2[s]';'glucose[s]'};
    expTask1.LBin=[23.6;23.6];
    expTask1.UBin=[23.8;23.8];
    expTask1.outputs={'H2O[s]';'CO2[s]'};
    expTask1.LBout=[26.1;26.1];
    expTask1.UBout=[26.2;26.2];
    expTask1.equations={'ATP[c] + H2O[c] => ADP[c] + Pi[c] + H+[c]'};
    expTask1.LBequ=30.2;
    expTask1.UBequ=30.6;
    expTask1.changed={'ATP[a] + H2O[a] => ADP[a] + Pi[a] + H+[a]'};
    expTask1.LBrxn=56.2;
    %check that the check works as expected
    verifyNotEqual(testCase,taskStruct(1),expTask1)
    expTask1.UBrxn=60;%now add the last

    
    verifyEqual(testCase,taskStruct(1),expTask1)
    verifyEqual(testCase,taskStructExcel(1),expTask1)
    
    %check that we have 2 tasks in total
    verifyEqual(testCase,length(taskStruct),2)
    verifyEqual(testCase,length(taskStructExcel),2)

    expTask2.id='BS';
    expTask2.description='ATP de novo synthesis';
    expTask2.shouldFail=false;
    expTask2.printFluxes=false;
    expTask2.comments='';
    expTask2.inputs={'O2[s]';'glucose[s]';'NH3[s]';'Pi[s]'};
    expTask2.LBin=[0;0;0;0];
    expTask2.UBin=[1000;1000;1000;1000];
    expTask2.outputs={'H2O[s]';'CO2[s]';'ATP[c]'};
    expTask2.LBout=[0;0;1];
    expTask2.UBout=[1.1;1.1;1.3];
    expTask2.equations={'ATP[c] + H2O[c] <=> ADP[c] + Pi[c] + H+[c]'};
    expTask2.LBequ=-1000;
    expTask2.UBequ=1000;
    expTask2.changed={};
    expTask2.LBrxn=[];
    expTask2.UBrxn=[];

    verifyEqual(testCase,taskStruct(2),expTask2)
    verifyEqual(testCase,taskStructExcel(2),expTask2)

    %and, check that the two formats produce the same
    verifyEqual(testCase,taskStruct,taskStructExcel)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModel = getTstModel()
    testModel = struct();
    testModel.id = 'testModel';
    testModel.rxns = {};
    testModel.S=[];
    testModel.rev=[];
    testModel.mets = {'a';'a2';'b';'c';'d';'e';'e2';'f'};
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
    rxnsToAdd.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};;
    testModel = addRxns(testModel,rxnsToAdd, 3);
    testModel.c = [0;0;0;0;0;0;0;1;0;0];%optimize for output flux, if this is used, not sure
    testModel.ub = repmat(1000,10,1);
    testModel.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
    testModel.rxnNames = testModel.rxns;
    testModel.b = repmat(0,8,1);
end

function testModelRxnScores = getTstModelRxnScores()
    testModelRxnScores = [-2;-2;-1;7;0.5;0.5;-1;-2;-1;1.5];
end

function testModelTasks = getTstModelTasks()
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
% testModel2 - not used directly for now, but indirectly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModel2 = getTstModel2()
    testModel2 = struct();
    testModel2.id = 'testModel2';
    testModel2.rxns = {};
    testModel2.S=[];
    testModel2.rev=[];
    testModel2.mets = {'a';'b'};
    testModel2.metNames = {'a';'b'};
    testModel2.comps = {'s'};
    testModel2.compNames = testModel2.comps;
    testModel2.metComps = [1;1];
    testModel2.genes = {'G1';'G2';'G3';'G4'};
    testModel2.grRules = {};
    testModel2.rxnGeneMat = [];

    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'R1';'R2';'R3';'R4'};
    rxnsToAdd.equations = {'a[s] <=>';...
                           'a[s] => b[s]';...
                           'a[s] <=> b[s]';...
                           'b[s] =>'};
    rxnsToAdd.grRules = testModel2.genes;
    testModel2 = addRxns(testModel2,rxnsToAdd, 3,true,true);
    testModel2.c = [0;0;0;1];%optimize for output flux, if this is used, not sure
    testModel2.ub = repmat(1000,4,1);
    testModel2.lb = [-1000;0;-1000;0];
    testModel2.rxnNames = testModel2.rxns;
    testModel2.b = zeros(2,1);
end

%    testRxnScores2 = [-1.1;-1;8;-1];
%    testRxnScores2b = [-1.1;8;-1;-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModel4 = getTstModel4()
    testModel2 = getTstModel2();
    testModel4 = testModel2;

    rxnsToAdd = struct();
    rxnsToAdd.rxns = {'R5';'R6';'R7';'R8';'R9';'R10';'R11'};
    rxnsToAdd.equations = {'5 a[s] <=> 5 d[s]';...
                           'e[s] <=> d[s]';
                           'f[s] + g[s] <=> e[s]';
                           'b[s] <=> f[s]';
                           'h[s] <=> g[s]';
                           'h[s] =>';
                           'e[s] => g[s]'};
    rxnsToAdd.grRules = {'G5';'G6';'G7';'G8';'G9';'G10';'G11'};
    testModel4 = addRxns(testModel4,rxnsToAdd, 3, [], true, true);
end

function testModel4RxnScores = getTstModel4RxnScores()
    testModel4RxnScores = [-1;-1;2;-1;0.5;-2;1;1.3;-0.5;-0.4;8];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0001: testModel without tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0001(testCase)
%    detectedMets = {};
    testParams = struct();

%    params.TimeLimit = 10;

    testModel = getTstModel();
    prepDataTest1 = prepINITModel(testModel, {}, {}, false);
    %check some things in the prepData
    %1. We expect 3 rxns in origRxnsToZero:
    verifyTrue(testCase, all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreImportRxns) , {'R1';'R2';'R8'})))
    %note that R7 should not be there, since it has a GPR.

    arrayData1.genes = testModel.genes;
    arrayData1.tissues = {'a'};
    arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
    arrayData1.threshold = 1;
    tst1ResModel1 = ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);

    %We expect R1, R2, and R8 to be added since they have no GPRs and are exch/simple transport rxns
    %R7 however will not be added since it has a GPR
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R8';'R9';'R10'})))

    %also test spontaneous
    prepDataTest1 = prepINITModel(testModel, {}, {'R7';'R10'}, false);
    verifyTrue(testCase, all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreImportRxns | prepDataTest1.toIgnoreSpont), {'R1';'R2';'R7';'R8';'R10'})))
    arrayData1.genes = testModel.genes;
    arrayData1.tissues = {'a'};
    arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
    arrayData1.threshold = 1;
    tst1ResModel1 = ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);
    %the model should now change to include the "correct" path and 
    %skip R9/R10:
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8'})), 1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0002: Create a task that wants to generate e[s] from a[s] for testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0002(testCase)
    testModel = getTstModel();
    testModelTasks = getTstModelTasks();
    testRxnScores = getTstModelRxnScores();
    testParams = struct();
    prepDataTest1 = prepINITModel(testModel, testModelTasks, {}, false);
    %We now expect to R2 and R7 to be essential. Note that R1 and R8 are not essential,
    %the exchange rxns are not used when checking tasks.
    %This is a bit complicated to check, because the essential rxns are expressed
    %as rxn ids of the linearly merged model. We expect those to be called R1 and R7:
    verifyTrue(testCase, all(strcmp(prepDataTest1.essentialRxns,{'R1';'R7'})))

    arrayData1.genes = testModel.genes;
    arrayData1.tissues = {'a'};
    arrayData1.levels = getExprForRxnScore(testRxnScores);
    arrayData1.threshold = 1;
    tst1ResModel1 = ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);
    %Since both R2 and R7 are now essential, we expect all rxns to be on except R3 and 
    %R5 (which have a negative total score and are not needed for the task)
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0003: The second step - gapfilling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0003(testCase)
    %First generate a model with gaps. We can use testModel. First we add 
    %boundary mets. Then we remove the exchange reactions 
    %(which is required for filling gaps) and create a gap by removing the R7 reaction.
    testModel = getTstModel();
    testModelTasks = getTstModelTasks();
    testRxnScores = getTstModelRxnScores();

    mTempRef = addBoundaryMets(testModel);
    mTempRef = removeReactions(mTempRef, {'R1';'R8'});
    mTemp = removeReactions(mTempRef, {'R7'});
    mTemp.id = 'tmp';
    tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
    %now check that R7 is added back
    [outModel,addedRxnMat] = fitTasksOpt(mTemp,mTempRef,[],true,min(tmpRxnScores,-0.1),testModelTasks);
    verifyTrue(testCase, all(strcmp(mTempRef.rxns(addedRxnMat),'R7')))%ok
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0004: MergeLinear and groupRxnScores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0004(testCase)
    %first testModel
    testModel = getTstModel();
    testRxnScores = getTstModelRxnScores();
    [reducedModel,origRxnIds,groupIds,reversedRxns]=mergeLinear(testModel, {});
    %We expect mergeLinear to merge {R1,R2}, {R3,R5}, {R4,R6}, {R7,R8}, {R9,R10}
    verifyTrue(testCase, all(groupIds == [1;1;2;3;2;3;4;4;5;5]))
    %we expect R1, R3, R4, R7 to be irreversible, R9 to be reversible
    verifyTrue(testCase, all(reducedModel.rev == [0;0;0;0;1]))
    verifyTrue(testCase, all(reducedModel.lb == [0;0;0;0;-1000]))

    newRxnScores=groupRxnScores(reducedModel, testRxnScores, origRxnIds, groupIds, ismember(origRxnIds, {'R1';'R2';'R8'}));
    verifyTrue(testCase, all(newRxnScores == [0;-0.5;7.5;-1;0.5]))

    %then testModel4
    testModel4 = getTstModel4();
    [reducedModel,origRxnIds,groupIds,reversedRxns]=mergeLinear(testModel4, {});
    %we expect {R5,R6},{R7,R8}, and{R9,R10} to be merged
    verifyTrue(testCase, all(groupIds == [0;0;0;0;1;1;2;2;3;3;0]))
    %check reversibility
    verifyTrue(testCase, all(reducedModel.rev == [1;0;1;0;1;1;0;0]))
    %check that some reactions have flipped direction when turned to irrev
    verifyTrue(testCase, strcmp(constructEquations(reducedModel, 'R9'),'g[s] => '))
    verifyTrue(testCase, all(find(reversedRxns) == [6;9]))
    %constructEquations(testModel4)
    %constructEquations(reducedModel)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0006: reverseRxns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testftINIT_T0006(testCase)
    %R1 = '=> a[s]';...
    %R3 = 'a[c] <=> b[c] + c[c]';...
    testModel = getTstModel();

    tmpModel = reverseRxns(testModel, {'R1';'R3'});
    res = constructEquations(tmpModel, {'R1';'R3'});
    expRes = {'a[s] => ';'b[c] + c[c] <=> a[c]'}
    verifyTrue(testCase, all(strcmp(res,expRes)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0007: rescaleModelForINIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0007(testCase)
    miniModel = struct();
    miniModel.S = [1,1000;-1,-1000];
    res = rescaleModelForINIT(miniModel,10);
    expRes = [1,10;-1,-10];
    verifyTrue(testCase, all(all(res.S == expRes))) %ok
end



