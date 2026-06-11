classdef tComparison < RavenTestCase
% tComparison  Tests for the multi-model comparison functions in comparison/.

    methods (Test)

        function compareRxnsGenesMetsCompsRuns(testCase)
            modelA = testCase.model; modelA.id = 'modelA';
            modelB = removeReactions(testCase.model, testCase.model.rxns(1:5), ...
                true, true, true); modelB.id = 'modelB';
            evalc('cs = compareRxnsGenesMetsComps({modelA, modelB});');
            testCase.verifyClass(cs, 'struct');
        end

        function compareMultipleModelsRuns(testCase)
            modelA = testCase.model; modelA.id = 'modelA';
            modelB = removeReactions(testCase.model, testCase.model.rxns(1:5), ...
                true, true, true); modelB.id = 'modelB';
            evalc('cs = compareMultipleModels({modelA, modelB});');
            testCase.verifyClass(cs, 'struct');
        end

    end
end
