classdef tScoringEquivalence < RavenTestCase
% tScoringEquivalence  scoreModel as a wrapper over scoreComplexModel.
%
%   scoreModel is the tINIT-era scorer and scoreComplexModel the ftINIT one.
%   They are now one implementation: scoreModel keeps its own argument order
%   and its -Inf convention for unmeasured genes, and delegates the scoring.
%
%   The Matches* tests are the contract the wrapper has to keep. The Changes*
%   tests pin the three behaviours the fold deliberately altered, so that none
%   of them can drift back unnoticed.

    methods (Test)

        %% What the wrapper preserves

        function arrayDataMaxScoringMatches(testCase)
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            simple  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1', ...
                'multipleGeneScoring', 'best', 'multipleCellScoring', 'best');
            complex = scoreComplexModel(m, [], a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'multipleCellScoring', 'max');
            testCase.verifyEqual(complex, simple, 'AbsTol', 1e-12);
        end

        function hpaDataMaxScoringMatches(testCase)
            m = testCase.scoringTestModel();
            h = testCase.scoringHpaData();
            simple  = scoreModel(m, h, 'tissue', 't1', ...
                'multipleGeneScoring', 'best', 'multipleCellScoring', 'best');
            complex = scoreComplexModel(m, h, [], 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'multipleCellScoring', 'max');
            testCase.verifyEqual(complex, simple, 'AbsTol', 1e-12);
        end

        function noGeneScoreIsHonouredByBoth(testCase)
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            noGene  = -7;
            simple  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1', ...
                'noGeneScore', noGene);
            complex = scoreComplexModel(m, [], a, 't1', 'noGeneScore', noGene, ...
                'isozymeScoring', 'max', 'complexScoring', 'max');
            noGeneRxn = strcmp(m.rxns, 'R1');
            testCase.verifyEqual(simple(noGeneRxn), noGene);
            testCase.verifyEqual(complex(noGeneRxn), noGene);
        end

        function unmeasuredGeneScoreStaysMinusInf(testCase)
            % The one convention the wrapper still translates: an unmeasured
            % gene is -Inf here and NaN in scoreComplexModel, because
            % removeLowScoreGenes reads NaN as "no evidence" and leaves the
            % gene alone, where -Inf < 0 would prune it.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            a.genes     = a.genes(1:2);        % drop G3 and G4
            a.levels    = a.levels(1:2, :);
            a.threshold = a.threshold(1:2);

            [~, gSimple]  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1');
            [~, gComplex] = scoreComplexModel(m, [], a, 't1');
            unmeasured = ~ismember(m.genes, a.genes);

            testCase.verifyTrue(all(isinf(gSimple(unmeasured))));
            testCase.verifyTrue(all(isnan(gComplex(unmeasured))));
            % Where a gene does have data the two agree exactly
            testCase.verifyEqual(gSimple(~unmeasured), gComplex(~unmeasured), ...
                'AbsTol', 1e-12);
        end

        function hpaArrayPrecedenceIsAnOption(testCase)
            % scoreModel decides per reaction: any gene of a reaction having
            % HPA data means array data is ignored for that whole reaction.
            % That is now dataPrecedence 'reaction' rather than a difference
            % between two implementations. R7 is "G1 or G4", G1 has HPA data
            % and G4 only array data.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            h = testCase.scoringHpaDataLowG1();
            simple = scoreModel(m, h, 'arrayData', a, 'tissue', 't1');
            perGene = scoreComplexModel(m, h, a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max');
            perRxn  = scoreComplexModel(m, h, a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'dataPrecedence', 'reaction');
            r7 = strcmp(m.rxns, 'R7');

            % Per reaction, only G1's HPA score is in scope
            testCase.verifyEqual(simple(r7), -8, 'AbsTol', 1e-12);
            testCase.verifyEqual(perRxn(r7), -8, 'AbsTol', 1e-12);
            % Per gene, G4's array score is seen too, and it is higher
            testCase.verifyEqual(perGene(r7), 5*log(3), 'AbsTol', 1e-12);
            % The wrapper picks the per-reaction rule, so it matches that one
            testCase.verifyEqual(simple, perRxn, 'AbsTol', 1e-12);
        end

        %% What the fold changed

        function changesAverageToReduceDownTheRule(testCase)
            % multipleGeneScoring 'average' used to be a flat mean over a
            % reaction's genes. It now averages down the grRule, so a complex
            % counts once against its isozymes rather than once per subunit.
            % R5 is "(G1 and G2) or G3". No caller uses 'average'.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            simple  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1', ...
                'multipleGeneScoring', 'average');
            complex = scoreComplexModel(m, [], a, 't1', ...
                'isozymeScoring', 'average', 'complexScoring', 'average');
            r5 = strcmp(m.rxns, 'R5');

            g = 5*log([1.5; 2; 2.5]);          % scores of G1, G2, G3
            testCase.verifyEqual(simple(r5), mean([mean(g(1:2)); g(3)]), ...
                'AbsTol', 1e-12);
            testCase.verifyNotEqual(simple(r5), mean(g));   % the old flat mean
            testCase.verifyEqual(simple(r5), complex(r5), 'AbsTol', 1e-12);
        end

        function changesMultiCellTypeReductionToSkipUnmeasured(testCase)
            % scoreModel used to reduce over a sparse genes-by-cell-types
            % matrix, so a cell type a gene was never measured in entered the
            % reduction as a numeric 0. Since 0 sits between 'Low' (10) and
            % 'Not detected' (-8), a gene not detected in one of two cell
            % types scored 0 under 'best' and was never pruned. The reduction
            % now runs over the measurements that exist.
            m = testCase.scoringTestModel();
            h = struct();
            h.genes      = {'G1'};
            h.tissues    = {'t1'; 't1'};
            h.celltypes  = {'ct1'; 'ct2'};
            h.levels     = {'High', 'Medium', 'Low', 'None'};
            h.gene2Level = [4 0];              % 'None' in ct1, unmeasured in ct2

            simple  = scoreModel(m, h, 'tissue', 't1', 'multipleCellScoring', 'best');
            complex = scoreComplexModel(m, h, [], 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'multipleCellScoring', 'max');
            r2 = strcmp(m.rxns, 'R2');         % grRule is just "G1"

            testCase.verifyEqual(simple(r2), -8, 'AbsTol', 1e-12);
            testCase.verifyEqual(complex(r2), -8, 'AbsTol', 1e-12);
        end

        function changesScoringSourceToGrRules(testCase)
            % scoreModel used to read rxnGeneMat and now reads grRules, like
            % scoreComplexModel. A model carrying its gene associations in
            % only one of the two now scores the same way through either
            % entry point. standardizeGrRules exists to keep them in step.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            noRules = m;
            noRules.grRules(:) = {''};         % rxnGeneMat left intact

            simple  = scoreModel(noRules, [], 'arrayData', a, 'tissue', 't1');
            complex = scoreComplexModel(noRules, [], a, 't1');
            r2 = strcmp(m.rxns, 'R2');

            testCase.verifyEqual(simple(r2), -2);        % the noGeneScore default
            testCase.verifyEqual(complex(r2), -2);
        end

        function rejectsUnknownDataPrecedence(testCase)
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            testCase.verifyError( ...
                @() scoreComplexModel(m, [], a, 't1', 'dataPrecedence', 'model'), ...
                'RAVEN:badInput');
        end

    end

    methods (Access = private)

        function m = scoringTestModel(~)
            % Small model whose grRules cover every rule shape the scoring
            % treats differently: none, single, OR, AND, nested, and a
            % mixed-data-source OR.
            m = struct();
            m.id = 'scoringTest';
            m.rxns = {}; m.S = []; m.rev = [];
            m.mets     = {'ac'; 'bc'; 'cc'; 'dc'; 'ec'};
            m.metNames = {'a'; 'b'; 'c'; 'd'; 'e'};
            m.comps = {'c'}; m.compNames = m.comps;
            m.metComps = [1; 1; 1; 1; 1];
            m.genes = {'G1'; 'G2'; 'G3'; 'G4'};
            m.grRules = {}; m.rxnGeneMat = [];
            r = struct();
            r.rxns = {'R1'; 'R2'; 'R3'; 'R4'; 'R5'; 'R6'; 'R7'};
            r.equations = {'=> a[c]'; 'a[c] => b[c]'; 'b[c] => c[c]'; ...
                'c[c] => d[c]'; 'a[c] => d[c]'; 'd[c] => e[c]'; 'e[c] =>'};
            r.grRules = {''; 'G1'; 'G1 or G2'; 'G1 and G2'; ...
                '(G1 and G2) or G3'; 'G4'; 'G1 or G4'};
            evalc('m = addRxns(m, r, 3);');
            m.c = zeros(7, 1);
            m.lb = zeros(7, 1);
            m.ub = repmat(1000, 7, 1);
            m.rxnNames = m.rxns;
            m.b = zeros(5, 1);
            evalc('[m.grRules, m.rxnGeneMat] = standardizeGrRules(m, true);');
        end

        function a = scoringArrayData(~)
            % Distinct, uncapped scores: with a threshold of 1 the score of
            % each gene is 5*log(level), i.e. 2.03, 3.47, 4.58 and 5.49.
            a = struct();
            a.genes     = {'G1'; 'G2'; 'G3'; 'G4'};
            a.tissues   = {'t1'; 't2'};
            a.celltypes = {'ct1'; 'ct2'};
            a.levels    = [1.5 1; 2 1; 2.5 1; 3 1];
            a.threshold = ones(4, 1);
        end

        function h = scoringHpaData(~)
            % One cell type, all four genes covered, one level each.
            h = struct();
            h.genes      = {'G1'; 'G2'; 'G3'; 'G4'};
            h.tissues    = {'t1'};
            h.celltypes  = {'ct1'};
            h.levels     = {'High', 'Medium', 'Low', 'None'};
            h.gene2Level = [1; 2; 3; 4];
        end

        function h = scoringHpaDataLowG1(~)
            % Only G1 has HPA data, and its level is the lowest one, so any
            % array score for the other genes outranks it.
            h = struct();
            h.genes      = {'G1'};
            h.tissues    = {'t1'};
            h.celltypes  = {'ct1'};
            h.levels     = {'High', 'Medium', 'Low', 'None'};
            h.gene2Level = 4;
        end

    end
end
