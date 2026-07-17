classdef tINIT < RavenTestCase
% tINIT  Tests for the ftINIT context-specific modelling functions in INIT/.
%
%   The MILP-based pipeline is exercised end-to-end on small synthetic models
%   (built by the local getTstModel* helpers, illustrated in tINITtestInfo/):
%   ftINITPipelineRuns / ftINITWithTaskRuns / ftINITMetabolomicsRuns /
%   ftINITFullVsThreeStepRuns cover prepINITModel + ftINIT and the internal
%   ftINITInternalAlg / groupRxnScores helpers, ftINITFillGapsForAllTasksRuns
%   covers the gap-filling step, all guarded on a MILP solver (Gurobi/SCIP).
%   The self-contained helpers (mergeLinear, groupRxnScores, reverseRxns,
%   rescaleModelForINIT, scoreComplexModel, getExprForRxnScore) are tested
%   directly against known results.

    methods (Test)

        function getINITStepsReturnsSteps(testCase)
            steps = getINITSteps([], '1+1');
            testCase.verifyNotEmpty(steps);
        end

        function getExprForRxnScoreRuns(testCase)
            expr = getExprForRxnScore(rand(10, 1), 1);
            testCase.verifyNotEmpty(expr);
        end

        function removeLowScoreGenesRuns(testCase)
            geneScores = randn(numel(testCase.model.genes), 1);
            evalc('m2 = removeLowScoreGenes(testCase.model, geneScores);');
            testCase.verifyClass(m2, 'struct');
        end

        function reverseRxnsRuns(testCase)
            % R1 = '=> a[s]'; R3 = 'a[c] <=> b[c] + c[c]'
            testModel = getTstModel();
            tmpModel  = reverseRxns(testModel, {'R1';'R3'});
            res       = constructEquations(tmpModel, {'R1';'R3'});
            expRes    = {'a[s] => ';'b[c] + c[c] <=> a[c]'};
            testCase.verifyTrue(all(strcmp(res, expRes)));
        end

        function rescaleModelForINITRuns(testCase)
            miniModel      = struct();
            miniModel.S    = [1,1000;-1,-40];
            miniModel.rxns = {'1';'2'};
            miniModel.mets = {'1';'2'};
            res = rescaleModelForINIT(miniModel, 10);
            testCase.verifyTrue(abs(res.S(1,2) - res.S(2,2)*-10) < 10^-6);
            testCase.verifyTrue(abs((abs(res.S(1,2)) + abs(res.S(2,2)))/2) - 1 < 10^-6);
        end

        function mergeLinearGroupsAndScores(testCase)
            % mergeLinear merges linear reaction chains; groupRxnScores then
            % aggregates the per-reaction scores onto the merged groups.
            testModel     = getTstModel();
            testRxnScores = getTstModelRxnScores();
            evalc('[reducedModel,origRxnIds,groupIds,reversedRxns] = mergeLinear(testModel, {});');
            % {R1,R2}, {R3,R5}, {R4,R6}, {R7,R8}, {R9,R10} are merged
            testCase.verifyTrue(all(groupIds == [1;1;2;3;2;3;4;4;5;5]));
            % R1, R3, R4, R7 irreversible, R9 reversible
            testCase.verifyTrue(all(reducedModel.rev == [0;0;0;0;1]));
            testCase.verifyTrue(all(reducedModel.lb == [0;0;0;0;-1000]));

            newRxnScores = groupRxnScores(reducedModel, testRxnScores, ...
                origRxnIds, groupIds, ismember(origRxnIds, {'R1';'R2';'R8'}));
            testCase.verifyTrue(all(newRxnScores == [0;-0.5;7.5;-1;0.5]));

            % testModel4 has reactions that are not all linearly merged
            testModel4 = getTstModel4();
            evalc('[reducedModel,origRxnIds,groupIds,reversedRxns] = mergeLinear(testModel4, {});');
            % {R5,R6}, {R7,R8}, {R9,R10} are merged
            testCase.verifyTrue(all(groupIds == [0;0;0;0;1;1;2;2;3;3;0]));
            testCase.verifyTrue(all(reducedModel.rev == [1;0;1;0;1;1;0;0]));
            % some reactions flip direction when turned irreversible
            testCase.verifyTrue(strcmp(constructEquations(reducedModel, 'R9'), 'g[s] => '));
            testCase.verifyTrue(all(find(reversedRxns) == [6;9]));
        end

        function scoreComplexModelRuns(testCase)
            arrayData.genes   = testCase.model.genes;
            arrayData.tissues = {'t1'};
            arrayData.levels  = abs(randn(numel(testCase.model.genes), 1)) + 1;
            arrayData.threshold = 1;
            evalc('rxnScores = scoreComplexModel(testCase.model, [], arrayData, ''t1'', []);');
            testCase.verifyNumElements(rxnScores, numel(testCase.model.rxns));
        end

        function scoreComplexModelHpaUnmeasuredCellTypeIsNotZero(testCase)
            % A gene not detected in one cell type of a tissue and not
            % measured in the other must keep its 'Not detected' score. An
            % unmeasured cell type is not a measurement of zero, and zero
            % ranks above 'Not detected' (-8), so scoring it as such would
            % also keep the gene from ever being pruned.
            m = testCase.model;
            hpaData.genes      = m.genes(1);
            hpaData.tissues    = {'liver';'liver'};
            hpaData.celltypes  = {'hepatocyte';'kupffer'};
            hpaData.levels     = {'Not detected'};
            hpaData.gene2Level = sparse(1,2);
            hpaData.gene2Level(1,1) = 1;   % not detected in hepatocyte
                                           % not measured in kupffer

            evalc(['[~,~,hpaScores] = scoreComplexModel(m, hpaData, [], ' ...
                '''liver'', ''multipleCellScoring'', ''max'');']);
            testCase.verifyEqual(hpaScores(1), -8, 'AbsTol', 1e-9);

            % The average is likewise taken over the measurements that exist
            evalc(['[~,~,hpaScores] = scoreComplexModel(m, hpaData, [], ' ...
                '''liver'', ''multipleCellScoring'', ''average'');']);
            testCase.verifyEqual(hpaScores(1), -8, 'AbsTol', 1e-9);
        end

        function scoreComplexModelExactScores(testCase)
            % scoreComplexModel must reproduce the known reaction scores for
            % the reference model prepared by prepINITModel.
            testCase.assumeMILPSolver();
            testModel = getTstModel();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;
            rxnScores = scoreComplexModel(prepData.refModel, [], arrayData, arrayData.tissues{1}, []);
            testCase.verifyTrue(all(abs(rxnScores - getTstModelRxnScores()) < 10^-10));
        end

        function ftINITPipelineRuns(testCase)
            % prepINITModel + ftINIT end-to-end on testModel without tasks.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            % R1 and R8 are exchange reactions without a GPR
            testCase.verifyTrue(all(strcmp( ...
                prepData.refModel.rxns(prepData.toIgnoreExch), {'R1';'R8'})));

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % R1 and R8 are added (no GPR, exchange); R2 (transport, no GPR) is
            % removed by the standard third step; R7 is not added (it has a GPR).
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R4';'R6';'R8';'R9';'R10'})));

            % make R7 and R10 spontaneous
            evalc('prepData = prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');
            testCase.verifyTrue(all(strcmp( ...
                prepData.refModel.rxns(prepData.toIgnoreExch | prepData.toIgnoreSpont), ...
                {'R1';'R7';'R8';'R10'})));
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % the model now takes the "correct" path (incl. R2), skipping R9/R10
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R4';'R6';'R7';'R8'})));
        end

        function ftINITWithTaskRuns(testCase)
            % A task requiring e[s] from a[s] forces R2 and R7 to be essential.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testParams     = struct();
            evalc('prepData = prepINITModel(testModel, testModelTasks, {}, false, {}, ''s'');');
            % essential rxns are reported on the linearly merged model (R1, R7)
            testCase.verifyTrue(all(strcmp(prepData.essentialRxns, {'R1';'R7'})));

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % R2 and R7 essential => all rxns on except R3 and R5 (negative
            % score, not needed for the task)
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
        end

        function ftINITFillGapsForAllTasksRuns(testCase)
            % Remove exchange rxns (required for gap filling) and create a gap
            % by removing R7; ftINITFillGapsForAllTasks must add R7 back.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R8'});
            mTemp    = removeReactions(mTempRef, {'R7'});
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
            evalc(['[~,addedRxnMat] = ftINITFillGapsForAllTasks(mTemp,mTempRef,' ...
                '[],false,min(tmpRxnScores,-0.1),testModelTasks);']);
            testCase.verifyTrue(all(strcmp(mTempRef.rxns(addedRxnMat), 'R7')));
        end

        function ftINITFillGapsForAllTasksReportsUnfillableTask(testCase)
            % R7 is required for the task, so removing it from the reference
            % model as well leaves the task unfillable. That must be reported,
            % not reported as "Added 0 reaction(s)" like an already-working
            % task.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R7';'R8'});
            mTemp    = mTempRef;
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;9;10]);
            evalc(['[~,addedRxnMat,failedTasks] = ftINITFillGapsForAllTasks(mTemp,' ...
                'mTempRef,[],false,min(tmpRxnScores,-0.1),testModelTasks);']);
            testCase.verifyTrue(failedTasks(1));
            testCase.verifyFalse(any(addedRxnMat(:)));
        end

        function ftINITMetabolomicsRuns(testCase)
            % Detected metabolites steer ftINIT towards alternative pathways.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            % adding metabolite f favours R9/R10 again (presence of R2 is random)
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''f''},getINITSteps(),true,true,testParams,false);']);
            if length(resModel.rxns) == 7
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R4';'R6';'R7';'R8';'R9';'R10'})));
            else
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
            end

            % metabolites a, e, f give the same result
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''f'';''a'';''e''},getINITSteps(),true,true,testParams,false);']);
            if length(resModel.rxns) == 7
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R4';'R6';'R7';'R8';'R9';'R10'})));
            else
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
            end

            % metabolite b turns on R2 and R3/R5 and turns off R9
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''b''},getINITSteps(),true,true,testParams,false);']);
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8'})));

            % testModel5 contains reactions that are not linearly merged
            testModel5 = getTstModel5();
            arrayData.genes  = testModel5.genes;
            arrayData.levels = getExprForRxnScore(getTstModel5RxnScores());

            evalc('prepData5 = prepINITModel(testModel5, {}, {''R7'';''R10''}, false, {}, ''s'');');
            evalc(['resModel = ftINIT(prepData5,arrayData.tissues{1},[],[],' ...
                'arrayData,{},getINITSteps(),true,true,testParams,false);']);
            % a->g->e via R11/R13 (score -2) ties the R2 path (score -2); the
            % solver takes the R11/R13 route, avoiding R9/R10.
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R4';'R6';'R7';'R8';'R11';'R13'})));

            % adding metabolite g drops R7
            evalc('prepData5 = prepINITModel(testModel5, {}, {''R10''}, false, {}, ''s'');');
            arrayData.levels(7) = getExprForRxnScore(-1.1); % avoid randomness
            evalc(['resModel = ftINIT(prepData5,arrayData.tissues{1},[],[],' ...
                'arrayData,{''g''},getINITSteps(),true,true,testParams,false);']);
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R4';'R6';'R8';'R11';'R13'})));
        end

        function ftINITFullVsThreeStepRuns(testCase)
            % On a larger model the three-step and "full" variants must give
            % consistent results.
            testCase.assumeMILPSolver();
            testModelL          = getTstModelL();
            testModelLGeneScores = getTstModelLGeneScores();
            testParams          = struct();

            arrayDataL.genes     = testModelL.genes;
            arrayDataL.tissues   = {'t1'};
            arrayDataL.levels    = getExprForRxnScore(testModelLGeneScores, 1);
            arrayDataL.threshold = 1;

            evalc('prepDataL = prepINITModel(testModelL, [], {}, false, {}, ''s'');');
            evalc(['mres = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],' ...
                'arrayDataL,[],getINITSteps(),true,true,testParams,false);']);
            evalc(['mres2 = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],' ...
                'arrayDataL,[],getINITSteps([], ''full''),true,true,testParams,false);']);

            expResult = {'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2'; ...
                'E3';'E4';'E5';'E6';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7'; ...
                'R8';'R9';'R12';'R13';'R14';'R15'};
            testCase.verifyTrue(all(contains(mres.rxns,  expResult)));
            testCase.verifyTrue(all(contains(mres2.rxns, expResult)));
        end

        function ftINITSeriesVariantsRun(testCase)
            % Only '1+1' and 'full' were ever exercised, so the 2-step series
            % from the paper -- and the allowExcretion constraint they lean on
            % -- had no coverage at all.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            steps = getINITSteps([], '2+1');
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],steps,true,true,testParams,false);']);
            % Same answer as the '1+1' series in ftINITPipelineRuns.
            testCase.verifyEqual(resModel.rxns, {'R1';'R4';'R6';'R8';'R9';'R10'});

            steps = getINITSteps([], '2+0');
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],steps,true,true,testParams,false);']);
            % '2+0' skips step 3, so the GPR-less transport R2 survives.
            testCase.verifyEqual(resModel.rxns, {'R1';'R2';'R4';'R6';'R8';'R9';'R10'});
        end

    end
end

%==========================================================================
% Synthetic test models, illustrated in tINITtestInfo/06.svg and 07.svg.
%==========================================================================

function testModel = getTstModel()
testModel = struct();
testModel.id = 'testModel';
testModel.rxns = {};
testModel.S = [];
testModel.rev = [];
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
evalc('testModel = addRxns(testModel,rxnsToAdd, 3);');
testModel.c = [0;0;0;0;0;0;0;1;0;0];
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

function testModel2 = getTstModel2()
testModel2 = struct();
testModel2.id = 'testModel2';
testModel2.rxns = {};
testModel2.S = [];
testModel2.rev = [];
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
evalc('testModel2 = addRxns(testModel2,rxnsToAdd,3,'''',true,true);');
testModel2.c = [0;0;0;1];
testModel2.ub = repmat(1000,4,1);
testModel2.lb = [-1000;0;-1000;0];
testModel2.rxnNames = testModel2.rxns;
testModel2.b = zeros(2,1);
end

function testModel4 = getTstModel4()
testModel4 = getTstModel2();

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

function testModelL = getTstModelL()
testModelL = struct();
testModelL.id = 'testModel';
testModelL.rxns = {};
testModelL.S = [];
testModelL.rev = [];
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
evalc('testModelL = addRxns(testModelL,rxnsToAdd, 3);');
testModelL.c = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
testModelL.rxnNames = testModelL.rxns;
testModelL.b = repmat(0,length(testModelL.mets),1);
end

function testModelLGeneScores = getTstModelLGeneScores()
%testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};
testModelLGeneScores = [3; -1;   8;    6;    -5;    5;     4;    5;    2;    3;    6;    1;    3;    1;    -3;    1;     3;      1;    2];
end
