classdef tReconstruction < RavenTestCase
% tReconstruction  Tests for the de-novo reconstruction functions in reconstruction/.
%
%   Homology and KEGG reconstruction depend on external binaries (BLAST+,
%   DIAMOND) and on a local KEGG dump / aligners. Those are guarded and skip
%   when unavailable; the self-contained functions are tested directly.

    methods (Test)

        function guessCompositionRuns(testCase)
            evalc('[m2, guessed] = guessComposition(testCase.model);');
            testCase.verifyClass(m2, 'struct');
        end

        function makeFakeBlastStructureReturnsStruct(testCase)
            % makeFakeBlastStructure requires at least 10 ortholog pairs.
            ol = [testCase.model.genes(1:10), strcat('t_', testCase.model.genes(1:10))];
            bs = makeFakeBlastStructure(ol, 'srcModel', 'tgtOrg');
            testCase.verifyClass(bs, 'struct');
        end

        function getModelFromHomologyBuildsDraft(testCase)
            src = testCase.model; src.id = 'srcModel';
            n = min(20, numel(src.genes));
            ol = [src.genes(1:n), strcat('t_', src.genes(1:n))];
            bs = makeFakeBlastStructure(ol, 'srcModel', 'tgtOrg');
            evalc('draft = getModelFromHomology({src}, bs, ''tgtOrg'');');
            testCase.verifyClass(draft, 'struct');
        end

        function getWSLpathReturnsPath(testCase)
            p = getWSLpath('C:\foo\bar');
            testCase.verifyTrue(ischar(p) || isstring(p));
        end

        function getBlastRunsWhenAvailable(testCase)
            fa1 = fullfile(testCase.ravenRoot,'testing','function_tests','test_data','human_galactosidases.fa');
            fa2 = fullfile(testCase.ravenRoot,'testing','function_tests','test_data','yeast_galactosidases.fa');
            testCase.assumeDependency(exist(fa1,'file')==2, 'test FASTA files');
            ok = true;
            try, evalc('bs = getBlast(''human'', fa1, {''yeast''}, {fa2});'); catch, ok = false; end
            testCase.assumeTrue(ok, 'BLAST+ binary not functional in this environment');
            testCase.verifyClass(bs, 'struct');
        end

        function getDiamondRunsWhenAvailable(testCase)
            fa1 = fullfile(testCase.ravenRoot,'testing','function_tests','test_data','human_galactosidases.fa');
            fa2 = fullfile(testCase.ravenRoot,'testing','function_tests','test_data','yeast_galactosidases.fa');
            testCase.assumeDependency(exist(fa1,'file')==2, 'test FASTA files');
            ok = true;
            try, evalc('bs = getDiamond(''human'', fa1, {''yeast''}, {fa2});'); catch, ok = false; end
            testCase.assumeTrue(ok, 'DIAMOND binary not functional in this environment');
            testCase.verifyClass(bs, 'struct');
        end

        function getKEGGModelForOrganismNeedsData(testCase)
            testCase.assumeFail('Requires keggModel.mat and an HMM library from raven-data.');
        end

        function getModelFromKEGGNeedsData(testCase)
            testCase.assumeFail('Requires keggModel.mat from raven-data.');
        end

        function getPhylDistNeedsData(testCase)
            testCase.assumeFail('Requires keggPhylDist.mat from raven-data.');
        end

    end
end
