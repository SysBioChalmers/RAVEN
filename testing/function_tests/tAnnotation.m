classdef tAnnotation < RavenTestCase
% tAnnotation  Tests for the annotation/metadata functions in annotation/.

    methods (Test)

        function editMiriamAddsAnnotation(testCase)
            m2 = editMiriam(testCase.model, 'met', 1, 'bigg.metabolite', 'testval', 'add');
            testCase.verifyClass(m2, 'struct');
            testCase.verifyEqual(numel(m2.mets), numel(testCase.model.mets));
        end

        function extractMiriamReturnsNames(testCase)
            [miriams, names] = extractMiriam(testCase.model.metMiriams); %#ok<ASGLU>
            testCase.verifyNotEmpty(names);
        end

        function assignSBOtermsRuns(testCase)
            evalc('m2 = assignSBOterms(testCase.model);');
            testCase.verifyClass(m2, 'struct');
        end

        function deltaGCsvRoundTrip(testCase)
            m = testCase.model;
            m.metDeltaG = (1:numel(m.mets))';
            m.rxnDeltaG = (1:numel(m.rxns))';
            metCsv = [tempname '.csv'];
            rxnCsv = [tempname '.csv'];
            testCase.addTeardown(@() delete(metCsv));
            testCase.addTeardown(@() delete(rxnCsv));
            evalc('deltaGCSV(m, ''save'', metCsv, rxnCsv);');
            base = testCase.model;
            evalc('m2 = deltaGCSV(base, ''load'', metCsv, rxnCsv);');
            testCase.verifyEqual(m2.metDeltaG, m.metDeltaG, 'AbsTol', 1e-6);
            testCase.verifyEqual(m2.rxnDeltaG, m.rxnDeltaG, 'AbsTol', 1e-6);
        end

    end
end
