classdef tComparison < RavenTestCase
% tComparison  Tests for the multi-model comparison functions in comparison/.

    methods (Test)

        function compareMultipleModelsRuns(testCase)
            modelA = testCase.model; modelA.id = 'modelA';
            modelB = removeReactions(testCase.model, testCase.model.rxns(1:5), ...
                true, true, true); modelB.id = 'modelB';
            evalc('cs = compareMultipleModels({modelA, modelB});');
            testCase.verifyClass(cs, 'struct');
        end

        function compareMultipleModelsReportsIdentityOverlap(testCase)
            % Covers the overlap comparison merged in from the now-removed
            % compareRxnsGenesMetsComps: rxns/mets/genes/eccodes/metNames/equ/uEqu.
            modelA = testCase.model; modelA.id = 'modelA';
            removedRxn = testCase.model.rxns{1};
            modelB = removeReactions(testCase.model, testCase.model.rxns(1), ...
                true, true, true); modelB.id = 'modelB';
            evalc('cs = compareMultipleModels({modelA, modelB});');
            for field = {'rxns','mets','genes','eccodes','metNames','equ','uEqu'}
                testCase.verifyTrue(isfield(cs, field{1}));
                testCase.verifyTrue(isfield(cs.(field{1}), 'comparison'));
                testCase.verifyTrue(isfield(cs.(field{1}), 'nElements'));
            end
            % modelA has one more reaction than modelB, so the rxns overlap must
            % show a combination where only modelA is included.
            onlyA = cs.rxns.comparison(:,1) & ~cs.rxns.comparison(:,2);
            testCase.verifyTrue(any(cs.rxns.nElements(onlyA) > 0));
            testCase.verifyTrue(any(strcmp(modelA.rxns, removedRxn)));
        end

        function diffModelsEqualToItself(testCase)
            r = diffModels(testCase.model, testCase.model);
            testCase.verifyTrue(r.equal);
            testCase.verifyEmpty(r.differences);
        end

        function diffModelsDetectsDroppedReaction(testCase)
            modelB = removeReactions(testCase.model, testCase.model.rxns(1), ...
                true, true, true);
            r = diffModels(testCase.model, modelB);
            testCase.verifyFalse(r.equal);
            testCase.verifyTrue(ismember(testCase.model.rxns{1}, r.rxnsOnlyInA));
            testCase.verifyTrue(any(contains(r.differences, 'reactions only in A')));
        end

        function diffModelsDetectsBoundDifference(testCase)
            modelB = testCase.model;
            modelB.lb(1) = modelB.lb(1) - 42;
            r = diffModels(testCase.model, modelB);
            testCase.verifyFalse(r.equal);
            testCase.verifyTrue(any(contains(r.differences, 'bounds')));
        end

        function diffModelsStoichWithinToleranceIsEqual(testCase)
            modelB = testCase.model;
            nz = find(modelB.S(:,1), 1);
            modelB.S(nz,1) = modelB.S(nz,1) + 1e-12;
            r = diffModels(testCase.model, modelB, 'stoichTol', 1e-9);
            testCase.verifyTrue(r.equal);
        end

        function diffModelsStoichOutsideToleranceIsNotEqual(testCase)
            modelB = testCase.model;
            nz = find(modelB.S(:,1), 1);
            modelB.S(nz,1) = modelB.S(nz,1) + 5;
            r = diffModels(testCase.model, modelB);
            testCase.verifyFalse(r.equal);
            testCase.verifyTrue(any(contains(r.differences, 'coef[')));
        end

        function diffModelsGrRuleIsOrderInsensitive(testCase)
            % The improvement over a string compare: "a and b" == "b and a".
            withGpr = find(cellfun(@(r) contains(r,' and '), testCase.model.grRules), 1);
            testCase.assumeNotEmpty(withGpr, 'model has an AND grRule');
            modelB = testCase.model;
            genes = strtrim(strsplit(modelB.grRules{withGpr}, ' and '));
            modelB.grRules{withGpr} = strjoin(fliplr(genes), ' and ');
            r = diffModels(testCase.model, modelB);
            % this reaction must not be reported as a grRule difference
            testCase.verifyFalse(any(contains(r.differences, ...
                [testCase.model.rxns{withGpr} ': grRule'])));
        end

        function diffModelsGrRuleLogicalChangeIsDetected(testCase)
            withGpr = find(cellfun(@(r) contains(r,' and '), testCase.model.grRules), 1);
            testCase.assumeNotEmpty(withGpr, 'model has an AND grRule');
            modelB = testCase.model;
            modelB.grRules{withGpr} = strrep(modelB.grRules{withGpr}, ' and ', ' or ');
            r = diffModels(testCase.model, modelB);
            testCase.verifyFalse(r.equal);
            testCase.verifyTrue(any(contains(r.differences, ...
                [testCase.model.rxns{withGpr} ': grRule'])));
        end

        function diffModelsTruncatesPerCategory(testCase)
            modelB = testCase.model;
            modelB.lb = modelB.lb - 1; % a bound diff on every reaction
            r = diffModels(testCase.model, modelB, 'maxPerCategory', 5);
            testCase.verifyFalse(r.equal);
            testCase.verifyTrue(any(contains(r.differences, 'truncated at 5')));
        end

    end
end
