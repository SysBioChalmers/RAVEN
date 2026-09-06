classdef tScoringEquivalence < RavenTestCase
% tScoringEquivalence  Can scoreComplexModel stand in for scoreModel?
%
%   scoreModel and scoreComplexModel both turn HPA and/or array data into a
%   score per reaction, but they are separate implementations: scoreModel
%   reduces over the genes in rxnGeneMat, scoreComplexModel evaluates the
%   grRule with a configurable operator for AND and for OR.
%
%   The tests below establish where the two agree and where they cannot be
%   made to agree by choice of options. The Matches* tests pin the settings
%   under which scoreComplexModel reproduces scoreModel; the Diverges* tests
%   pin the behaviours that a replacement would have to account for, and are
%   expected to fail if either function is changed to close the gap.
%
%   Agreement requires isozymeScoring = complexScoring = 'max', which makes
%   every grRule collapse to the maximum over its genes, i.e. what
%   scoreModel computes with multipleGeneScoring = 'best'.

    methods (Test)

        %% Where the two agree

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
            % A reaction without genes takes noGeneScore in both functions,
            % and the option is respected rather than hard-coded.
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

        function geneScoresMatchWhereDataExists(testCase)
            % The second output is built the same way in both functions for
            % every gene that actually has a measurement.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            [~, gSimple]  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1');
            [~, gComplex] = scoreComplexModel(m, [], a, 't1');
            measured = ismember(m.genes, a.genes);
            testCase.verifyEqual(gComplex(measured), gSimple(measured), 'AbsTol', 1e-12);
        end

        %% Where they cannot be made to agree

        function hpaArrayPrecedenceDiverges(testCase)
            % scoreModel decides per reaction: if any gene of a reaction has
            % HPA data, array data is ignored for that whole reaction.
            % scoreComplexModel decides per gene: array scores are written
            % first and HPA overwrites only the genes it covers. R7 is
            % "G1 or G4", G1 has HPA data and G4 only array data, so the two
            % see a different set of scores for the same rule.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            h = testCase.scoringHpaDataLowG1();
            simple  = scoreModel(m, h, 'arrayData', a, 'tissue', 't1', ...
                'multipleGeneScoring', 'best', 'multipleCellScoring', 'best');
            complex = scoreComplexModel(m, h, a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'multipleCellScoring', 'max');
            r7 = strcmp(m.rxns, 'R7');

            % scoreModel sees only G1's HPA score
            testCase.verifyEqual(simple(r7), -8, 'AbsTol', 1e-12);
            % scoreComplexModel also sees G4's array score, which is higher
            testCase.verifyEqual(complex(r7), 5*log(3), 'AbsTol', 1e-12);
            testCase.verifyNotEqual(complex(r7), simple(r7));
        end

        function averageScoringDivergesOnNestedRules(testCase)
            % 'average' means a flat mean over the reaction's genes in
            % scoreModel, but a mean of means down the rule tree in
            % scoreComplexModel. R5 is "(G1 and G2) or G3".
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            simple  = scoreModel(m, [], 'arrayData', a, 'tissue', 't1', ...
                'multipleGeneScoring', 'average');
            complex = scoreComplexModel(m, [], a, 't1', ...
                'isozymeScoring', 'average', 'complexScoring', 'average');
            r5 = strcmp(m.rxns, 'R5');

            g = 5*log([1.5; 2; 2.5]);          % scores of G1, G2, G3
            testCase.verifyEqual(simple(r5), mean(g), 'AbsTol', 1e-12);
            testCase.verifyEqual(complex(r5), mean([mean(g(1:2)); g(3)]), 'AbsTol', 1e-12);
            testCase.verifyNotEqual(complex(r5), simple(r5));
        end

        function multiCellTypeHpaReductionDiverges(testCase)
            % When a gene is measured in some but not all cell types,
            % scoreModel reduces over a sparse genes-by-cell-types matrix,
            % so unmeasured pairs enter the reduction as a numeric 0. A
            % negative level therefore reduces to 0 under 'best'.
            % scoreComplexModel reduces only over the cell types the gene
            % was actually measured in.
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

            testCase.verifyEqual(simple(r2), 0, 'AbsTol', 1e-12);
            testCase.verifyEqual(complex(r2), -8, 'AbsTol', 1e-12);
            testCase.verifyNotEqual(complex(r2), simple(r2));
        end

        function scoringSourceFieldDiverges(testCase)
            % scoreModel reads rxnGeneMat, scoreComplexModel reads grRules.
            % A model carrying gene associations in only one of the two is
            % scored by one function and falls through to noGeneScore in the
            % other.
            m = testCase.scoringTestModel();
            a = testCase.scoringArrayData();
            noRules = m;
            noRules.grRules(:) = {''};         % rxnGeneMat left intact

            simple  = scoreModel(noRules, [], 'arrayData', a, 'tissue', 't1');
            complex = scoreComplexModel(noRules, [], a, 't1');
            r2 = strcmp(m.rxns, 'R2');

            testCase.verifyEqual(simple(r2), 5*log(1.5), 'AbsTol', 1e-12);
            testCase.verifyEqual(complex(r2), -2);       % the noGeneScore default
            testCase.verifyNotEqual(complex(r2), simple(r2));
        end

        function unmeasuredGeneScoreDiverges(testCase)
            % The geneScores output marks a gene without data as -Inf in
            % scoreModel and as NaN in scoreComplexModel. Callers that test
            % for one will not see the other.
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
        end

    end

    methods (Access = private)

        function m = scoringTestModel(~)
            % Small model whose grRules cover every rule shape the two
            % scoring functions treat differently: none, single, OR, AND,
            % nested, and a mixed-data-source OR.
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
