%run this test case with the command
%results = runtests('tinitTests.m')
function tests = tinitTests
tests = functiontests(localfunctions);
solverExist = [exist('gurobi','file'), exist('scip.mexw64','file')] ==3;
if solverExist(1)
    try
        gurobi_read('solverTests.m'); % Random function call
    catch ME
        if ~startsWith(ME.message,'Gurobi error 10012') % Expected error code, others may indicate problems with license
            solverExist(1) = false;
        end
    end
end
if all(~solverExist)
    disp('No suitable solver (gurobi or scip) was found, ftINIT tests skipped.')
    skipTests = contains({tests.Name},'ftinit','IgnoreCase',true);
    tests(skipTests) = [];
end
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

function testModelRxnScores = getTstModelRxnScores()
testModelRxnScores = [-2;-2;-1;7;0.5;0.5;-1;-2;-3;3.5];
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
testModel2 = addRxns(testModel2,rxnsToAdd,3,'',true,true);
testModel2.c = [0;0;0;1];%optimize for output flux, if this is used, not sure
testModel2.ub = repmat(1000,4,1);
testModel2.lb = [-1000;0;-1000;0];
testModel2.rxnNames = testModel2.rxns;
testModel2.b = zeros(2,1);
end

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
[~,testModel4] = evalc("addRxns(testModel4,rxnsToAdd, 3, [], true, true);");
end

function testModel4RxnScores = getTstModel4RxnScores()
testModel4RxnScores = [-1;-1;2;-1;0.5;-2;1;1.3;-0.5;-0.4;8];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testModel5 = getTstModel5()
testModel = getTstModel();

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R11';'R12';'R13';'R14'};
rxnsToAdd.equations = {'a[c] <=> g[c]';...
    'a[c] <=> g[c]';...
    'g[c] <=> e[c]';...
    'g[c] <=> e[c]'};
rxnsToAdd.grRules = {'G11';'G12';'G13';'G14'};
[~,testModel5] = evalc("addRxns(testModel,rxnsToAdd, 3, [], true, true);");
end

function testModel5RxnScores = getTstModel5RxnScores()
testModel5RxnScores = [getTstModelRxnScores();-1;-1.5;-1;-1.5];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0001: testModel without tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0001(testCase)
solverExist = [exist('gurobi','file'), exist('scip.mexw64','file')] == 3;
currSolver = getpref('RAVEN','solver');
if solverExist(1)
    setRavenSolver('gurobi');
elseif solverExist(2)
    setRavenSolver('scip');
else
    error('No compatible solvers found');
end
%    detectedMets = {};
testParams = struct();

%    params.TimeLimit = 10;

testModel = getTstModel();
[~, prepDataTest1] = evalc('prepINITModel(testModel, {}, {}, false, {}, ''s'');');
%check some things in the prepData
%1. We expect 3 rxns in origRxnsToZero:
verifyTrue(testCase, all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch) , {'R1';'R8'})))
%note that R7 should not be there, since it has a GPR.

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
arrayData1.threshold = 1;

[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],getINITSteps(),true,true,testParams,false);');

%We expect R1 and R8 to be added since they have no GPRs and are exch rxns. The transport rxn R2 without GPR will however be removed,
%since we in the standard setting run the third step with [1;0;0;0;1;0;0], meaning that such reactions will be removed
%R7 however will not be added since it has a GPR
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R4';'R6';'R8';'R9';'R10'})))

%also test spontaneous
[~, prepDataTest1] = evalc('prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');
verifyTrue(testCase, all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreSpont), {'R1';'R7';'R8';'R10'})))
arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
arrayData1.threshold = 1;
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],getINITSteps(),true,true,testParams,false);');
%the model should now change to include the "correct" path (including 'R2') and
%skip R9/R10:
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8'})), 1)
setRavenSolver(currSolver);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0002: Create a task that wants to generate e[s] from a[s] for testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0002(testCase)
solverExist = [exist('gurobi','file'), exist('scip.mexw64','file')] ==3;
currSolver = getpref('RAVEN','solver');
if solverExist(1)
    setRavenSolver('gurobi');
elseif solverExist(2)
    setRavenSolver('scip');
else
    error('No compatible solvers found');
end
testModel = getTstModel();
testModelTasks = getTstModelTasks();
testRxnScores = getTstModelRxnScores();
testParams = struct();
[~, prepDataTest1] = evalc('prepINITModel(testModel, testModelTasks, {}, false, {}, ''s'');');
%We now expect to R2 and R7 to be essential. Note that R1 and R8 are not essential,
%the exchange rxns are not used when checking tasks.
%This is a bit complicated to check, because the essential rxns are expressed
%as rxn ids of the linearly merged model. We expect those to be called R1 and R7:
verifyTrue(testCase, all(strcmp(prepDataTest1.essentialRxns,{'R1';'R7'})))

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(testRxnScores);
arrayData1.threshold = 1;
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],getINITSteps(),true,true,testParams,false);');
%Since both R2 and R7 are now essential, we expect all rxns to be on except R3 and
%R5 (which have a negative total score and are not needed for the task)
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})))
setRavenSolver(currSolver);
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

mTempRef = closeModel(testModel);
mTempRef = removeReactions(mTempRef, {'R1';'R8'});
mTemp = removeReactions(mTempRef, {'R7'});
mTemp.id = 'tmp';
tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
%now check that R7 is added back
[~,outModel,addedRxnMat] = evalc('ftINITFillGapsForAllTasks(mTemp,mTempRef,[],false,min(tmpRxnScores,-0.1),testModelTasks);');
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
expRes = {'a[s] => ';'b[c] + c[c] <=> a[c]'};
verifyTrue(testCase, all(strcmp(res,expRes)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0007: rescaleModelForINIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0007(testCase)
miniModel = struct();
miniModel.S = [1,1000;-1,-40];
miniModel.rxns = {'1';'2'};
miniModel.mets = {'1';'2'};
res = rescaleModelForINIT(miniModel,10);
verifyTrue(testCase, abs(res.S(1,2) - res.S(2,2)*-10) < 10^-6)
verifyTrue(testCase, abs((abs(res.S(1,2)) + abs(res.S(2,2)))/2) - 1 < 10^-6)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0008: testModel with metabolomics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0008(testCase)
solverExist = [exist('gurobi','file'), exist('scip.mexw64','file')] ==3;
currSolver = getpref('RAVEN','solver');
if solverExist(1)
    setRavenSolver('gurobi');
elseif solverExist(2)
    setRavenSolver('scip');
else
    error('No compatible solvers found');
end

testParams = struct();

testModel = getTstModel();
[~, prepDataTest1] = evalc('prepINITModel(testModel, {}, {}, false, {}, ''s'');');

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
arrayData1.threshold = 1;

[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],getINITSteps(),true,true,testParams,false);');

%First the same as in T0001 (we keep them here to make the test case understandable):
%We expect R1 and R8 to be added since they have no GPRs and are exch rxns. The transport rxn R2 without GPR will however be removed,
%since we in the standard setting run the third step with [1;0;0;0;1;0;0], meaning that such reactions will be removed
%R7 however will not be added since it has a GPR
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R4';'R6';'R8';'R9';'R10'})))

%make R7 and R10 spontaneous (also same as in T0001)
[~, prepDataTest1] = evalc('prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');
verifyTrue(testCase, all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreSpont), {'R1';'R7';'R8';'R10'})))
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],getINITSteps(),true,true,testParams,false);');
%the model should now change to include the "correct" path (including 'R2') and
%skip R9/R10:
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8'})))
%now, test to add the metabolite f - this should turn back the favor to R9/R10:
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,{''f''},getINITSteps(),true,true,testParams,false);');
%R9 should now be included, the presence of R2 is random
if length(tst1ResModel1.rxns) == 7
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R4';'R6';'R7';'R8';'R9';'R10'})))
else
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})))
end
%now, test to add the metabolite a, e, f - this should give the same result:
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,{''f'';''a'';''e''},getINITSteps(),true,true,testParams,false);');
%Should be the same as above (R2 is random)
if length(tst1ResModel1.rxns) == 7
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R4';'R6';'R7';'R8';'R9';'R10'})))
else
    verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})))
end

%now, test to add the metabolite b - this should turn on R2 and R3/R5 and turn off R9:
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,{''b''},getINITSteps(),true,true,testParams,false);');
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8'})))

%now on model 5 to test reactions that are not linearly merged (testModel has only merged rxns)
testModel5 = getTstModel5();
arrayData1.genes = testModel5.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(getTstModel5RxnScores());
arrayData1.threshold = 1;

[~, prepDataTest5] = evalc('prepINITModel(testModel5, {}, {''R7'';''R10''}, false, {}, ''s'');');
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest5,arrayData1.tissues{1},[],[],arrayData1,{},getINITSteps(),true,true,testParams,false);');
%We expect the 'true' path, i.e. through R2, not R9/R10 or R11-R14
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8'})), 1)

%now add metabolite g
%modify the scores a bit
[~, prepDataTest5] = evalc('prepINITModel(testModel5, {}, {''R10''}, false, {}, ''s'');');
arrayData1.levels(7) = getExprForRxnScore(-1.1); %modify to avoid randomness
[~,tst1ResModel1] = evalc('ftINIT(prepDataTest5,arrayData1.tissues{1},[],[],arrayData1,{''g''},getINITSteps(),true,true,testParams,false);');
%We expect R2 to be replaced with R11 and R13
verifyTrue(testCase, all(strcmp(tst1ResModel1.rxns,{'R1';'R4';'R6';'R8';'R11';'R13'})), 1)
setRavenSolver(currSolver);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0009: getExprFromRxnScore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0009(testCase)
testModel = getTstModel();
[~, prepDataTest1] = evalc('prepINITModel(testModel, {}, {}, false, {}, ''s'');');

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(getTstModelRxnScores());
arrayData1.threshold = 1;

rxnScores = scoreComplexModel(prepDataTest1.refModel,[],arrayData1,arrayData1.tissues{1},[]);
expRes = getTstModelRxnScores();
verifyTrue(testCase, all(abs(rxnScores - expRes) < 10^-10)) %ok
end

function testModelL = getTstModelL()
testModelL = struct();
testModelL.id = 'testModel';
testModelL.rxns = {};
testModelL.S=[];
testModelL.rev=[];
testModelL.metNames = {'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'x1';'x2';'x3';'x4';'x5';'x6';'x7';'x8';'x9';'x10';'x11'};
testModelL.comps = {'s';'c'};
testModelL.compNames = testModelL.comps;
testModelL.metComps = [1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];
testModelL.mets = strcat(testModelL.metNames, testModelL.comps(testModelL.metComps));

testModelL.grRules = {};
testModelL.rxnGeneMat = [];

testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};

testModelL.ub = [];
testModelL.lb = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {  'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2';'E2b';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10';'R11';'R12';'R13';'R14';'R15'};
rxnsToAdd.grRules = {'';  '';  '';  '';  '';  '';  '';  '';  ''; 'Ge1';'Ge2';'';'';'Ge4';'Ge5';'';'Ge7';'';'Ge9'; 'Gr1';'Gr2';'Gr3';'';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'';'Gr14';'Gr15'};
rxnsToAdd.equations = {'e1[s] <=>';...
    'e2[s] <=>';...
    'e3[s] <=>';...
    'e4[s] <=>';...
    'e5[s] <=>';...
    'e6[s] <=>';...
    'e7[s] <=>';...
    'e8[s] <=>';...
    'e9[s] <=>';...
    'e1[s] <=> e1[c]';...
    'e2[s] <=> e2[c]';...
    'e2[s] <=> e2[c]';... %b variant
    'e3[s] <=> e3[c]';...
    'e4[s] <=> e4[c]';...
    'e5[s] <=> e5[c]';...
    'e6[s] <=> e6[c]';...
    'e7[s] <=> e7[c]';...
    'e8[s] <=> e8[c]';...
    'e9[s] <=> e9[c]';...
    'e1[c] + e2[c] <=> x1[c]';... %R1
    'e1[c] + e3[c] => x2[c] + x3[c]';... %R2
    'e4[c] + x3[c] => x4[c] + x5[c]';... %R3
    'e5[c] + e6[c] + x4[c] => 2 x2[c] + x6[c]';... %R4
    'x1[c] + x2[c] <=> x7[c] + 2 x8[c]';... %R5
    'x2[c] + x8[c] => x3[c] + x9[c]';... %R6
    'x4[c] <=> x9[c]';... %R7
    'x5[c] <=> x9[c]';... %R8
    'x6[c] <=> x10[c]';... %R9
    'x6[c] <=> x11[c]';... %R10
    'x10[c] + 2 x11[c] => e7[c]';... %R11
    'x9[c] + x10[c] <=> e8[c]';... %R12
    'x7[c] + x8[c] + x9[c] => e9[c]';... %R13
    'x6[c] => x9[c]';... %R14
    'x3[c] => x9[c]'... %R15
    };
testModelL = addRxns(testModelL,rxnsToAdd, 3);
testModelL.c = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];%optimize for output flux, if this is used, not sure
testModelL.rxnNames = testModelL.rxns;
testModelL.b = repmat(0,length(testModelL.mets),1);
end

function testModelLGeneScores = getTstModelLGeneScores()
%testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};
testModelLGeneScores = [3; -1;   8;    6;    -5;    5;     4;    5;    2;    3;    6;    1;    3;    1;    -3;    1;     3;      1;    2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0050: Test ftINIT on a slightly larger and more complex model
%       Specifically tests if the three-step variant and the "full"
%       variant gives similar results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testftINIT_T0050(testCase)
solverExist = [exist('gurobi','file'), exist('scip.mexw64','file')] ==3;
currSolver = getpref('RAVEN','solver');
if solverExist(1)
    setRavenSolver('gurobi');
elseif solverExist(2)
    setRavenSolver('scip');
else
    error('No compatible solvers found');
end

testModelL = getTstModelL();
testModelLGeneScores = getTstModelLGeneScores();
testParams = struct();

arrayDataL = struct();
arrayDataL.genes = testModelL.genes;
arrayDataL.tissues = {'t1'};
arrayDataL.levels = getExprForRxnScore(testModelLGeneScores,1);
arrayDataL.threshold = 1;

%Run prep data
[~, prepDataL] = evalc('prepINITModel(testModelL, [], {}, false, {}, ''s'');');

[~,mres] = evalc('ftINIT(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],getINITSteps(),true,true,testParams,false);');
[~,mres2] = evalc('ftINIT(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],getINITSteps([], ''full''),true,true,testParams,false);');

expResult = {  'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2';'E3';'E4';'E5';'E6';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R12';'R13';'R14';'R15'};

verifyTrue(testCase, all(contains(mres.rxns,expResult)))
verifyTrue(testCase, all(contains(mres2.rxns,expResult)))

%run the old tINIT version (in Human-GEM)
%this is just to show that they become different, not really part of the test case
%paramsL2 = struct();
%paramsL2.TimeLimit = 1000;
%testModelL2 = closeModel(testModelL);
%init_modelOrig = getINITModel2(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);

%in this call, I have modified the code - the possibility to turn off met secretion + don't allow flux in both directions is not possible.
%The following line, around line 337, is changed
%from:
%[~, deletedRxnsInINIT, metProduction] = runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,true,false,params);
%to:
%[~, deletedRxnsInINIT, metProduction] = runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,false,true,params);
%init_modelOrigNoSecrOneDirOnly = getINITModel2(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);

%The models init_modelOrigNoSecrOneDirOnly and mres2 are very similar, (only one exch rxn differ, which is expected)
%init_modelOrig is quite different, with a lot of gaps, and worse. So, the conclusion is that the new version does a pretty good job.
setRavenSolver(currSolver);
end
