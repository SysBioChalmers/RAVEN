classdef tConversion < RavenTestCase
% tConversion  Tests for the model format-conversion functions in conversion/.

    methods (Test)

        function identifierPrefixRoundTrip(testCase)
            % Adding then removing SBML identifier prefixes is the identity.
            [mp, ~] = addIdentifierPrefix(testCase.model);
            [mr, ~] = removeIdentifierPrefix(mp);
            testCase.verifyEqual(mr.rxns, testCase.model.rxns);
            testCase.verifyEqual(mr.mets, testCase.model.mets);
            testCase.verifyEqual(mr.genes, testCase.model.genes);
        end

        function ravenCobraWrapperMarksCobra(testCase)
            cobra = ravenCobraWrapper(testCase.model);
            testCase.verifyTrue(isfield(cobra, 'rules'));   % COBRA-only field
        end

        function ravenCobraWrapperRoundTrip(testCase)
            cobra = ravenCobraWrapper(testCase.model);
            back = ravenCobraWrapper(cobra);
            testCase.verifyEqual(back.rxns, testCase.model.rxns);
            testCase.verifyNumElements(back.mets, numel(testCase.model.mets));
            testCase.verifyEqual(size(back.S), size(testCase.model.S));
        end

        function standardizeFieldOrderPreservesFields(testCase)
            m2 = standardizeModelFieldOrder(testCase.model);
            testCase.verifyEqual(sort(fieldnames(m2)), sort(fieldnames(testCase.model)));
        end

        function standardizeFieldOrderIsIdempotent(testCase)
            m2 = standardizeModelFieldOrder(testCase.model);
            m3 = standardizeModelFieldOrder(m2);
            testCase.verifyEqual(fieldnames(m3), fieldnames(m2));
        end

    end
end
