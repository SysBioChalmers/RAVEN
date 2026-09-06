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

        function removeIdentifierPrefixSkipsMissingFields(testCase)
            % A model missing an optional field (genes, metNames, rxnNames,
            % id) must not make the default field list throw; it should
            % just be skipped.
            m = rmfield(testCase.model, {'genes','metNames','rxnNames','id'});
            m.rxns = strcat('R_', m.rxns);
            m2 = removeIdentifierPrefix(m);
            testCase.verifyFalse(any(startsWith(m2.rxns, 'R_')));
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

        function ravenCobraWrapperGeneFieldsUseGeneMiriams(testCase)
            % Gene annotation fields must come from geneMiriams, not from
            % whatever extractMiriam last returned for metMiriams.
            m = testCase.model;
            m.metMiriams = cell(numel(m.mets),1);
            m.metMiriams{1}.name  = {'hmdb'};
            m.metMiriams{1}.value = {'HMDB00001'};
            m.geneMiriams = cell(numel(m.genes),1);
            m.geneMiriams{1}.name  = {'ncbigene'};
            m.geneMiriams{1}.value = {'12345'};
            evalc('newModel = ravenCobraWrapper(m);');
            testCase.verifyEqual(numel(newModel.geneEntrezID), numel(newModel.genes));
            testCase.verifyEqual(newModel.geneEntrezID{1}, '12345');
        end

        function ravenCobraWrapperEscapesRegexGeneNames(testCase)
            % A gene id containing a regex metacharacter ('.') must be
            % matched literally in grRules, not as a wildcard that could
            % also match an unrelated gene id (e.g. 'G.1' vs 'GX1').
            m = testCase.model;
            m.genes{1} = 'G.1';
            m.genes{2} = 'GX1';
            m.grRules{3} = 'GX1';
            evalc('cobra = ravenCobraWrapper(m);');
            testCase.verifyEqual(cobra.rules{3}, 'x(2)');
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
