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

        function guessCompositionAccountsForStoichiometricCoefficient(testCase)
            % A -> 2 B, with A's formula known (CH4). B's own coefficient is
            % 2, so its atom counts must be half of A's (C0.5H4/2 = C0.5H2),
            % not identical to A's.
            model = tReconstruction.toyModel();
            model.metFormulas = {'CH4'; ''};
            evalc('[m2, guessedFor, couldNotGuess] = guessComposition(model);');
            testCase.verifyEqual(guessedFor, {'B'});
            testCase.verifyEmpty(couldNotGuess);
            testCase.verifyEqual(m2.metFormulas{2}, 'C0.5H2');
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

        function getModelFromHomologyComplexPolicyAndOptions(testCase)
            % A small template: a complex, an isozyme pair, a gene-free
            % reaction, and two template genes that map to one new gene.
            t = struct();
            t.id = 'tmpl'; t.name = 'tmpl';
            t.mets = {'a_c';'b_c';'d_c'}; t.metNames = {'A';'B';'D'};
            t.comps = {'c'}; t.compNames = {'cytosol'}; t.metComps = [1;1;1];
            t.rxns = {'Rcplx';'Riso';'Rfree';'Rdup'};
            t.rxnNames = t.rxns;
            t.S = sparse([-1 0 0 -1; 0 -1 0 1; 1 1 -1 0]);
            t.lb = [0;0;0;0]; t.ub = [1000;1000;1000;1000]; t.c = [0;0;0;0]; t.b = [0;0;0];
            t.rev = [0;0;0;0];
            t.grRules = {'tg1 and tg2'; 'tg3 or tg4'; ''; 'tg3 or tg4'};
            t.genes = {'tg1';'tg2';'tg3';'tg4'};
            t.rxnGeneMat = sparse([1 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 1 1]);
            t.rxnNotes = {'n1';'n2';'n3';'n4'};
            t.rxnConfidenceScores = [4;4;4;4];
            t.version = '9.9.9';
            % makeFakeBlastStructure wants at least 10 pairs; the extra ones
            % are genes the template does not have
            pairs = [{'tg1','ng1'; 'tg3','ng3'; 'tg4','ng3'}; ...
                     [strcat('x', string(1:7))', strcat('nx', string(1:7))']];
            pairs = cellstr(pairs);
            bs = makeFakeBlastStructure(pairs, 'tmpl', 'new');
            evalc('flag = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0);');
            evalc('keep = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0, ''complexPolicy'', ''keep'', ''keepGeneFree'', true, ''preserveNotes'', true);');
            evalc('drop = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0, ''complexPolicy'', ''drop'');');
            % flag: unmapped subunit kept as OLD_; gene-free reaction dropped
            testCase.verifyTrue(contains(flag.grRules{strcmp(flag.rxns,'Rcplx')}, 'OLD_tmpl_tg2'));
            testCase.verifyFalse(ismember('Rfree', flag.rxns));
            % two template genes mapping to one new gene appear once
            testCase.verifyEqual(flag.grRules{strcmp(flag.rxns,'Riso')}, 'ng3');
            % keep: subunit dropped, gene-free kept, template notes kept
            testCase.verifyEqual(keep.grRules{strcmp(keep.rxns,'Rcplx')}, 'ng1');
            testCase.verifyTrue(ismember('Rfree', keep.rxns));
            testCase.verifyEqual(keep.rxnNotes{strcmp(keep.rxns,'Rcplx')}, 'n1');
            testCase.verifyEqual(keep.rxnConfidenceScores(strcmp(keep.rxns,'Rcplx')), 4);
            % drop: the incomplete complex is removed
            testCase.verifyFalse(ismember('Rcplx', drop.rxns));
            testCase.verifyTrue(ismember('Riso', drop.rxns));
            % the template version is not inherited
            testCase.verifyFalse(isfield(flag, 'version') && strcmp(flag.version, '9.9.9'));
            testCase.verifyEqual(flag.rxnNotes{1}, 'Included by getModelFromHomology');
        end

        function getModelFromHomologyRejectsUnknownComplexPolicy(testCase)
            src = testCase.model; src.id = 'srcModel';
            bs = makeFakeBlastStructure([src.genes(1:10), strcat('t_', src.genes(1:10))], 'srcModel', 'tgtOrg');
            testCase.verifyError(@() getModelFromHomology({src}, bs, 'tgtOrg', 'complexPolicy', 'maybe'), ?MException);
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
            matFile = fullfile(testCase.ravenRoot, 'reconstruction', 'kegg', 'keggModel.mat');
            testCase.assumeFalse(exist(matFile, 'file') == 2, ...
                'keggModel.mat is present; skipping the build-from-artefacts path.');
            testCase.assumeFail('Downloads and assembles the full KEGG artefact set from raven-data; not run automatically.');
        end

        function buildGlobalGPRJoinsGenesThroughSharedKO(testCase)
            % Offline unit test of the join at the core of buildGlobalKEGGModel:
            % two reactions sharing a KO, a KO used by two organisms' genes, and
            % a reaction with no KO at all (should end up with no genes).
            rxns = {'R1'; 'R2'; 'R3'};
            koReaction = table({'K1'; 'K1'; 'K2'}, {'R1'; 'R2'; 'R2'}, ...
                'VariableNames', {'ko', 'reaction'});
            organismGeneKO = table({'a'; 'a'; 'b'}, {'g1'; 'g2'; 'g1'}, {'K1'; 'K2'; 'K1'}, ...
                'VariableNames', {'organism', 'gene', 'ko'});
            [genes, rxnGeneMat] = buildGlobalGPR(rxns, koReaction, organismGeneKO);
            testCase.verifyEqual(genes, {'a:g1'; 'a:g2'; 'b:g1'});
            expected = [1 0 1; 1 1 1; 0 0 0];
            testCase.verifyEqual(full(rxnGeneMat), expected);
        end

        function getPhylDistNeedsData(testCase)
            distFile = fullfile(testCase.ravenRoot, 'reconstruction', 'kegg', 'keggPhylDist.mat');
            testCase.assumeFalse(exist(distFile, 'file') == 2, ...
                'keggPhylDist.mat is present; skipping error-path test.');
            testCase.verifyError(@() getPhylDist(), 'getPhylDist:noData');
        end

    end

    methods (Static)
        function model = toyModel()
            % Single reaction: A -> 2 B, one compartment.
            model = struct();
            model.id          = 'toy';
            model.name        = 'toy';
            model.mets        = {'m1'; 'm2'};
            model.metNames    = {'A'; 'B'};
            model.comps       = {'c'};
            model.compNames   = {'cytosol'};
            model.metComps    = [1; 1];
            model.metFormulas = {''; ''};
            model.rxns        = {'R1'};
            model.rxnNames    = {'R1'};
            model.S           = sparse([-1; 2]);
            model.lb          = 0;
            model.ub          = 1000;
            model.rev         = 0;
            model.c           = 0;
            model.b           = zeros(2, 1);
        end
    end
end
