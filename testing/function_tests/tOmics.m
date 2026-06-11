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

        function parseHPAparsesSyntheticDump(testCase)
            % parseHPA reads a tab-separated HPA dump with the columns
            % Gene / Gene name / Tissue / Cell type / Level / Reliability.
            % Synthesise a tiny one so the parser is exercised without the
            % multi-gigabyte external download.
            g = testCase.model.genes(1:3);
            f = [tempname '.tsv'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            fprintf(fid, 'Gene\tGene name\tTissue\tCell type\tLevel\tReliability\n');
            fprintf(fid, '%s\t%s\tliver\thepatocyte\tHigh\tApproved\n', g{1}, g{1});
            fprintf(fid, '%s\t%s\tliver\thepatocyte\tMedium\tApproved\n', g{2}, g{2});
            fprintf(fid, '%s\t%s\tbrain\tneuron\tLow\tSupported\n', g{3}, g{3});
            fclose(fid);
            evalc('hpa = parseHPA(f);');
            testCase.verifyClass(hpa, 'struct');
            testCase.verifyEqual(sort(hpa.genes), sort(g));
            testCase.verifyNumElements(hpa.tissues, 2);   % liver+brain cell types
        end

        function parseHPArnaParsesSyntheticDump(testCase)
            % parseHPArna (version 19) reads Gene / Gene name / Tissue /
            % TPM / pTPM / NX and builds a genes-by-tissues TPM matrix.
            g = testCase.model.genes(1:3);
            f = [tempname '.tsv'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            fprintf(fid, 'Gene\tGene name\tTissue\tTPM\tpTPM\tNX\n');
            fprintf(fid, '%s\t%s\tliver\t12.5\t10.1\t8.0\n', g{1}, g{1});
            fprintf(fid, '%s\t%s\tbrain\t3.2\t2.9\t1.5\n', g{2}, g{2});
            fprintf(fid, '%s\t%s\tliver\t0.0\t0.0\t0.0\n', g{3}, g{3});
            fclose(fid);
            evalc('rna = parseHPArna(f);');
            testCase.verifyClass(rna, 'struct');
            testCase.verifyNumElements(rna.genes, 3);
            testCase.verifyEqual(size(rna.levels), [3 2]);   % 3 genes x 2 tissues
        end

    end
end
