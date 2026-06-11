classdef tOmics < RavenTestCase
% tOmics  Tests for the omics-integration functions in omics/.

    methods (Test)

        function scoreModelFromArrayData(testCase)
            arrayData.genes = testCase.model.genes;
            arrayData.tissues = {'tissue1'};
            arrayData.levels = abs(randn(numel(testCase.model.genes), 1)) + 1;
            arrayData.threshold = 1;
            evalc('rxnScores = scoreModel(testCase.model, [], arrayData, ''tissue1'', []);');
            testCase.verifyNumElements(rxnScores, numel(testCase.model.rxns));
        end

        function parseHPArequiresFile(testCase)
            f = fullfile(testCase.ravenRoot,'testing','unit_tests','test_data','hpa.xml');
            testCase.assumeDependency(exist(f,'file')==2, 'HPA xml dump');
            evalc('hpa = parseHPA(f);');
            testCase.verifyClass(hpa, 'struct');
        end

        function parseHPArnaRequiresFile(testCase)
            f = fullfile(testCase.ravenRoot,'testing','unit_tests','test_data','hpa_rna.tsv');
            testCase.assumeDependency(exist(f,'file')==2, 'HPA RNA dump');
            evalc('rna = parseHPArna(f);');
            testCase.verifyClass(rna, 'struct');
        end

    end
end
