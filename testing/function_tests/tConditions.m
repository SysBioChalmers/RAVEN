classdef tConditions < RavenTestCase
% tConditions  Tests for the condition-application functions in conditions/.

    methods (Test)

        function applyConditionSetsBounds(testCase)
            rxn = testCase.model.rxns{1};
            condition.bounds = {struct('rxn', rxn, 'lb', 0, 'ub', 0)};
            evalc('m2 = applyCondition(testCase.model, condition);');
            idx = strcmp(m2.rxns, rxn);
            testCase.verifyEqual(m2.ub(idx), 0);
            testCase.verifyEqual(m2.lb(idx), 0);
        end

    end
end
