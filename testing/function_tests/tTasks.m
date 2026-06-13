classdef tTasks < RavenTestCase
% tTasks  Tests for the metabolic-task functions in tasks/.

    methods (Test)

        function checkTasksSucceedsOnFullModel(testCase)
            m = testCase.taskTestModel();
            task = testCase.taskTestStruct();
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 1);
        end

        function checkTasksFailsWithoutKeyReaction(testCase)
            m = removeReactions(testCase.taskTestModel(), {'R2'});
            task = testCase.taskTestStruct();
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 0);
        end

        function parseTaskListReturnsStruct(testCase)
            taskFile = fullfile(testCase.ravenRoot,'testing','function_tests', ...
                'test_data','test_tasks.txt');
            testCase.assumeDependency(exist(taskFile,'file')==2, 'test_tasks.txt');
            ts = parseTaskList(taskFile);
            testCase.verifyClass(ts, 'struct');
            testCase.verifyNotEmpty(ts);
        end

        function parseTaskListMatchesExpected(testCase)
            % The text and Excel task definitions must parse to identical
            % structs with exactly the expected fields and values.
            dataDir   = fullfile(testCase.ravenRoot,'testing','function_tests','test_data');
            txtFile   = fullfile(dataDir,'test_tasks.txt');
            xlsFile   = fullfile(dataDir,'test_tasks.xls');
            testCase.assumeDependency(exist(txtFile,'file')==2, 'test_tasks.txt');
            testCase.assumeDependency(exist(xlsFile,'file')==2, 'test_tasks.xls');

            taskStruct      = parseTaskList(txtFile);
            taskStructExcel = parseTaskList(xlsFile);

            expTask1.id          = 'ER';
            expTask1.description = 'Aerobic rephosphorylation of ATP from glucose';
            expTask1.shouldFail  = true;
            expTask1.printFluxes = true;
            expTask1.comments    = 'Messed up reaction';
            expTask1.inputs      = {'O2[s]';'glucose[s]'};
            expTask1.LBin        = [23.6;23.6];
            expTask1.UBin        = [23.8;23.8];
            expTask1.outputs     = {'H2O[s]';'CO2[s]'};
            expTask1.LBout       = [26.1;26.1];
            expTask1.UBout       = [26.2;26.2];
            expTask1.equations   = {'ATP[c] + H2O[c] => ADP[c] + Pi[c] + H+[c]'};
            expTask1.LBequ       = 30.2;
            expTask1.UBequ       = 30.6;
            expTask1.changed     = {'ATP[a] + H2O[a] => ADP[a] + Pi[a] + H+[a]'};
            expTask1.LBrxn       = 56.2;
            % incomplete struct must not match
            testCase.verifyNotEqual(taskStruct(1), expTask1);
            expTask1.UBrxn = 60;
            testCase.verifyEqual(taskStruct(1), expTask1);
            testCase.verifyEqual(taskStructExcel(1), expTask1);

            testCase.verifyEqual(length(taskStruct), 2);
            testCase.verifyEqual(length(taskStructExcel), 2);

            expTask2.id          = 'BS';
            expTask2.description = 'ATP de novo synthesis';
            expTask2.shouldFail  = false;
            expTask2.printFluxes = false;
            expTask2.comments    = '';
            expTask2.inputs      = {'O2[s]';'glucose[s]';'NH3[s]';'Pi[s]'};
            expTask2.LBin        = [0;0;0;0];
            expTask2.UBin        = [1000;1000;1000;1000];
            expTask2.outputs     = {'H2O[s]';'CO2[s]';'ATP[c]'};
            expTask2.LBout       = [0;0;1];
            expTask2.UBout       = [1.1;1.1;1.3];
            expTask2.equations   = {'ATP[c] + H2O[c] <=> ADP[c] + Pi[c] + H+[c]'};
            expTask2.LBequ       = -1000;
            expTask2.UBequ       = 1000;
            expTask2.changed     = {};
            expTask2.LBrxn       = [];
            expTask2.UBrxn       = [];
            testCase.verifyEqual(taskStruct(2), expTask2);
            testCase.verifyEqual(taskStructExcel(2), expTask2);

            % the two formats must produce identical results
            testCase.verifyEqual(taskStruct, taskStructExcel);
        end

    end
end
