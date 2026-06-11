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
            taskFile = fullfile(testCase.ravenRoot,'testing','unit_tests', ...
                'test_data','test_tasks.txt');
            testCase.assumeDependency(exist(taskFile,'file')==2, 'test_tasks.txt');
            ts = parseTaskList(taskFile);
            testCase.verifyClass(ts, 'struct');
            testCase.verifyNotEmpty(ts);
        end

    end
end
