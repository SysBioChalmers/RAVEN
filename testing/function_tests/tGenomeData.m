classdef tGenomeData < RavenTestCase
% tGenomeData  Tests for the genome-data utilities used in homology-based
%              reconstruction (downloadGenomeData, getGeneData,
%              processProteinFastaFile, renameModelGenes).
%
%   Parsing and renaming are tested offline with small GFF3/FASTA fixtures
%   (one eukaryote, one prokaryote). The NCBI download is guarded and skips
%   when there is no internet connection.

    methods (Access = protected)
        function d = fixtureDir(testCase)
            d = fullfile(testCase.ravenRoot,'testing','unit_tests', ...
                'test_data','genomeData');
        end
    end

    methods (Test)

        function getGeneDataReadsEukaryote(testCase)
            % gene -> mRNA -> CDS: the protein accession must come from the
            % CDS protein_id (NP_), not the mRNA transcript (NM_), and the
            % two CDS lines of one protein collapse to a single row.
            T = getGeneData(fullfile(testCase.fixtureDir(),'euk.gff'));
            testCase.verifyEqual(height(T), 2);
            i1 = find(strcmp(T.locus_tag,'LT1'),1);
            testCase.verifyEqual(T.old_locus_tag{i1}, 'OLD1');
            testCase.verifyEqual(T.GeneID{i1}, '111');
            testCase.verifyEqual(T.gene_name{i1}, 'AAA');
            testCase.verifyEqual(T.GenBank_protein{i1}, 'NP_111.1');
            i2 = find(strcmp(T.locus_tag,'LT2'),1);
            testCase.verifyEqual(T.GenBank_protein{i2}, 'NP_222.1');
        end

        function getGeneDataReadsProkaryote(testCase)
            % gene -> CDS directly.
            T = getGeneData(fullfile(testCase.fixtureDir(),'prok.gff'));
            testCase.verifyEqual(height(T), 1);
            testCase.verifyEqual(T.locus_tag{1}, 'b0001');
            testCase.verifyEqual(T.GeneID{1}, '333');
            testCase.verifyEqual(T.gene_name{1}, 'ccc');
            testCase.verifyEqual(T.GenBank_protein{1}, 'WP_333.1');
        end

        function getGeneDataWritesNoFileByDefault(testCase)
            prevFiles = dir('geneData.tsv');
            getGeneData(fullfile(testCase.fixtureDir(),'prok.gff'));
            testCase.verifyEqual(numel(dir('geneData.tsv')), numel(prevFiles));
        end

        function processProteinFastaFileRenamesHeaders(testCase)
            T   = getGeneData(fullfile(testCase.fixtureDir(),'euk.gff'));
            out = processProteinFastaFile(fullfile(testCase.fixtureDir(), ...
                'euk_protein.faa'), T, 'locus_tag', tempdir);
            testCase.addTeardown(@() delete(out));
            fa = readFasta(out);
            headers = {fa.Header};
            % NP_111.1 -> LT1, NP_222.1 -> LT2, NP_999.9 has no match and is kept
            testCase.verifyTrue(ismember('LT1', headers));
            testCase.verifyTrue(ismember('LT2', headers));
            testCase.verifyTrue(any(startsWith(headers, 'NP_999.9')));
        end

        function renameModelGenesRenamesGenesAndGrRules(testCase)
            model.id         = 'test';
            model.rxns       = {'r1'; 'r2'};
            model.genes      = {'LT1'; 'LT2'};
            model.grRules    = {'LT1'; 'LT1 or LT2'};
            model.rxnGeneMat = sparse([1 0; 1 1]);

            T = getGeneData(fullfile(testCase.fixtureDir(),'euk.gff'));
            evalc('m2 = renameModelGenes(model, T, ''locus_tag'', ''gene_name'');');

            testCase.verifyEqual(sort(m2.genes), {'AAA'; 'BBB'});
            testCase.verifyFalse(any(contains(m2.grRules, 'LT1')));
            testCase.verifyTrue(any(contains(m2.grRules, 'AAA')));
            testCase.verifyTrue(any(contains(m2.grRules, 'BBB')));
            testCase.verifyEqual(size(m2.rxnGeneMat,2), numel(m2.genes));
        end

        function downloadGenomeDataDownloads(testCase)
            % Network test: download a small genome from NCBI Datasets.
            % Skipped when there is no internet connection.
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d,'s'));
            try
                [gff, faa] = downloadGenomeData('GCF_000005845.2', d, false);
            catch ME
                testCase.assumeFail(['Skipped, likely no internet connection: ' ME.message]);
                return
            end
            testCase.verifyTrue(isfile(gff));
            testCase.verifyTrue(isfile(faa));
        end

    end
end
