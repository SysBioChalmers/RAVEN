classdef tAnnotation < RavenTestCase
% tAnnotation  Tests for the annotation/metadata functions in annotation/.

    methods (Test)

        function editMiriamAddsAnnotation(testCase)
            m2 = editMiriam(testCase.model, 'met', 1, 'bigg.metabolite', 'testval', 'add');
            testCase.verifyClass(m2, 'struct');
            testCase.verifyEqual(numel(m2.mets), numel(testCase.model.mets));
        end

        function editMiriamRemovesAnnotation(testCase)
            m2 = editMiriam(testCase.model, 'met', 1, 'bigg.metabolite', 'testval', 'add');
            m3 = editMiriam(m2, 'met', 1, 'bigg.metabolite', '', 'remove');
            testCase.verifyClass(m3, 'struct');
            [~, names] = extractMiriam(m3.metMiriams(1));
            testCase.verifyFalse(any(strcmp('bigg.metabolite', names)));
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

        function deltaGCsvDropsSentinelByDefault(testCase)
            % Default missingValue is 1e7, matching raven_toolbox's
            % DELTA_G_MISSING: yeast-GEM's "no measurement" sentinel loads as
            % NaN, not as a real deltaG.
            m = testCase.model;
            metCsv = [tempname '.csv'];
            testCase.addTeardown(@() delete(metCsv));
            fid = fopen(metCsv, 'w');
            fprintf(fid, 'Var1,Var2\n%s,-42.5\n%s,10000000.0\n', m.mets{1}, m.mets{2});
            fclose(fid);
            evalc('md = deltaGCSV(m, ''load'', ''metCsv'', metCsv);');
            testCase.verifyEqual(md.metDeltaG(1), -42.5, 'AbsTol', 1e-9);
            testCase.verifyTrue(isnan(md.metDeltaG(2)));
        end

        function deltaGCsvKeepsSentinelWhenDisabled(testCase)
            % Passing [] disables the sentinel (raven_toolbox's missing_value
            % =None), storing the 1e7 value literally.
            m = testCase.model;
            metCsv = [tempname '.csv'];
            testCase.addTeardown(@() delete(metCsv));
            fid = fopen(metCsv, 'w');
            fprintf(fid, 'Var1,Var2\n%s,10000000.0\n', m.mets{1});
            fclose(fid);
            evalc('mk = deltaGCSV(m, ''load'', ''metCsv'', metCsv, ''missingValue'', []);');
            testCase.verifyEqual(mk.metDeltaG(1), 1e7, 'AbsTol', 1e-3);
        end

    end
end
