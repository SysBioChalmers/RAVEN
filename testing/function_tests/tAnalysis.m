classdef tAnalysis < RavenTestCase
% tAnalysis  Tests for the flux-analysis and simulation functions in analysis/.

    methods (Test)

        function getAllowedBoundsSizes(testCase)
            [mn, mx] = getAllowedBounds(testCase.model, 1:10);
            testCase.verifyNumElements(mn, 10);
            testCase.verifyNumElements(mx, 10);
            testCase.verifyTrue(all(mx >= mn - 1e-9));
        end

        function getEssentialRxnsReturnsList(testCase)
            evalc('[ess, idx] = getEssentialRxns(testCase.model);');
            testCase.verifyClass(ess, 'cell');
            testCase.verifyTrue(all(ismember(ess, testCase.model.rxns)));
        end

        function haveFluxReturnsLogical(testCase)
            I = haveFlux(testCase.model);
            testCase.verifyNumElements(I, numel(testCase.model.rxns));
        end

        function getMinNrFluxesReturnsFlux(testCase)
            testCase.assumeMILPSolver();
            evalc('[x, I, exitFlag] = getMinNrFluxes(testCase.model, testCase.model.rxns);');
            testCase.verifyNotEmpty(x);
        end

        function getAllSubGraphsReturnsResult(testCase)
            sg = getAllSubGraphs(testCase.model);
            testCase.verifyNotEmpty(sg);
        end

        function findGeneDeletionsRuns(testCase)
            evalc('[genes, fluxes] = findGeneDeletions(testCase.model, ''sgd'', ''fba'');');
            testCase.verifyNotEmpty(genes);
        end

        function followFluxesRuns(testCase)
            sol = solveLP(testCase.model);
            out = evalc('followFluxes(testCase.model, sol.x, 0, 1000);');
            testCase.verifyClass(out, 'char');
        end

        function followChangedRuns(testCase)
            % Aerobic vs anaerobic guarantees several reactions change, which
            % avoids the empty-selection edge case.
            solA = solveLP(testCase.model);
            o2exch = find(strcmp(testCase.model.rxnNames, 'O2 exchange'));
            modelAna = setParam(testCase.model, 'eq', testCase.model.rxns(o2exch), 0);
            solB = solveLP(modelAna);
            out = evalc('followChanged(testCase.model, solA.x, solB.x, 10, 0.01, 0);');
            testCase.verifyClass(out, 'char');
        end

        function getFluxZComputesScores(testCase)
            n = numel(testCase.model.rxns);
            Z = getFluxZ(rand(n, 20), rand(n, 20));
            testCase.verifyNumElements(Z, n);
        end

        function analyzeSamplingRuns(testCase)
            n = numel(testCase.model.rxns);
            Tex = randn(n, 1);           % per-reaction expression t-scores
            solutionsA = rand(n, 20);
            solutionsB = rand(n, 20);
            evalc('scores = analyzeSampling(Tex, 18, solutionsA, solutionsB, false);');
            testCase.verifySize(scores, [n 3]);
        end

        function randomSamplingRuns(testCase)
            evalc('sols = randomSampling(testCase.model, 5);');
            testCase.verifyEqual(size(sols, 1), numel(testCase.model.rxns));
        end

        function randomSamplingDeterministicWithSeed(testCase)
            % Two runs with the same seed must produce identical solution matrices
            % (SAMP2: seed parameter).  runParallel=false keeps the RNG state
            % sequential and avoids parallel-worker stream divergence.
            evalc('s1 = randomSampling(testCase.model, 5, ''seed'', 42, ''runParallel'', false);');
            evalc('s2 = randomSampling(testCase.model, 5, ''seed'', 42, ''runParallel'', false);');
            testCase.verifyEqual(s1, s2);
        end

        function reporterMetabolitesRuns(testCase)
            pvals = rand(numel(testCase.model.genes), 1);
            rm = reporterMetabolites(testCase.model, testCase.model.genes, pvals);
            testCase.verifyClass(rm, 'struct');
        end

        function reporterMetabolitesIsDeterministic(testCase)
            % Closed-form background correction (RM1) must produce identical
            % Z-scores on repeated calls with identical inputs.
            pvals = rand(numel(testCase.model.genes), 1);
            rm1 = reporterMetabolites(testCase.model, testCase.model.genes, pvals);
            rm2 = reporterMetabolites(testCase.model, testCase.model.genes, pvals);
            testCase.verifyEqual(rm1.metZScores, rm2.metZScores);
        end

        function FSEOFRuns(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            [~, exchIdx] = getExchangeRxns(testCase.model);
            targetRxn = testCase.model.rxns{exchIdx(1)};
            outFile = [tempname '.txt'];
            evalc('targets = FSEOF(testCase.model, biomassRxn, targetRxn, 5, 0.9, outFile);');
            testCase.verifyClass(targets, 'struct');
        end

        function runRobustnessAnalysisRuns(testCase)
            [~, exchIdx] = getExchangeRxns(testCase.model);
            controlRxn = testCase.model.rxns{exchIdx(1)};
            evalc('[cFlux, oFlux] = runRobustnessAnalysis(testCase.model, controlRxn, 5);');
            testCase.verifyNumElements(cFlux, 5);
        end

        function runProductionEnvelopeRuns(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            [~, exchIdx] = getExchangeRxns(testCase.model);
            targetRxn = testCase.model.rxns{exchIdx(1)};
            evalc('[bio, tgt] = runProductionEnvelope(testCase.model, targetRxn, biomassRxn, 5);');
            testCase.verifyNotEmpty(bio);
        end

        function runPhenotypePhasePlaneRuns(testCase)
            [~, exchIdx] = getExchangeRxns(testCase.model);
            r1 = testCase.model.rxns{exchIdx(1)};
            r2 = testCase.model.rxns{exchIdx(2)};
            evalc('g = runPhenotypePhasePlane(testCase.model, r1, r2, 3, 3);');
            testCase.verifyNotEmpty(g);
        end

        function runDynamicFBARuns(testCase)
            % Use glucose (consumed) as substrate and plot target; the function
            % always plots, so the plotted reaction must carry non-zero flux.
            [~, exchIdx] = getExchangeRxns(testCase.model);
            glcMet = find(strcmp(testCase.model.metNames, 'D-Glucose'), 1);
            glcExch = exchIdx(any(testCase.model.S(glcMet, exchIdx) ~= 0, 1));
            subRxn = testCase.model.rxns(glcExch);
            evalc(['c = runDynamicFBA(testCase.model, subRxn, 10, 0.1, 0.1, 2, ' ...
                'subRxn, {});']);
            testCase.verifyNotEmpty(c);
        end

        function runSimpleOptKnockRuns(testCase)
            testCase.assumeMILPSolver();
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            [~, exchIdx] = getExchangeRxns(testCase.model);
            targetRxn = testCase.model.rxns{exchIdx(1)};
            evalc(['out = runSimpleOptKnock(testCase.model, targetRxn, biomassRxn, ' ...
                'testCase.model.rxns(1:5), ''rxns'', 1);']);
            testCase.verifyClass(out, 'struct');
        end

    end
end
