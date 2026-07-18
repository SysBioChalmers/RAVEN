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

        function cycleFreeFluxRemovesLoop(testCase)
            % R1/R2 form a futile a<->b loop; the input carries 1000 around it
            % with a net through-flux of 10. cycleFreeFlux must strip the loop
            % (R1 -> 10, R2 -> 0) while preserving the exchange fluxes.
            m = tAnalysis.loopModel();
            cf = cycleFreeFlux(m, [10; 1000; 990; 10]);
            testCase.verifyEqual(cf(1), 10, 'AbsTol', 1e-6);   % R_in preserved
            testCase.verifyEqual(cf(4), 10, 'AbsTol', 1e-6);   % R_out preserved
            testCase.verifyEqual(cf(2), 10, 'AbsTol', 1e-6);   % R1 net only
            testCase.verifyLessThan(abs(cf(3)), 1e-6);         % R2 loop flux gone
        end

        function looplessFVAExcludesLoopFlux(testCase)
            % Standard FVA lets R1 reach 1000 around the loop; loopless FVA
            % must cap it at the net 10, and R2 (only ever loop flux) at 0.
            testCase.assumeMILPSolver();
            m = tAnalysis.loopModel();
            [lo, hi] = looplessFVA(m, {'R1';'R2'});
            testCase.verifyEqual(hi(1), 10, 'AbsTol', 1e-4);
            testCase.verifyLessThan(abs(hi(2)), 1e-4);
            testCase.verifyGreaterThanOrEqual(lo(1), -1e-4);
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

        function traceFluxPathReturnsCells(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            % Use two reactions that must be metabolically connected
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            [~, exchIdx] = getExchangeRxns(testCase.model);
            % Pick an exchange reaction with non-zero flux
            exchFlux = sol.x(exchIdx);
            active   = exchIdx(abs(exchFlux) > 1e-8);
            if isempty(active), testCase.assumeTrue(false); end
            targetRxn = testCase.model.rxns{active(1)};
            evalc('[p, m, f] = traceFluxPath(testCase.model, sol.x, biomassRxn, targetRxn, ''verbose'', false);');
            % Either a path is found (cell + double) or not — both are valid
            testCase.verifyClass(f, 'double');
            testCase.verifyGreaterThanOrEqual(f, 0);
            testCase.verifyLessThanOrEqual(f, 1);
        end

        function traceFluxPathSameReactionIsIdentity(testCase)
            testCase.assumeSolver('solveLP');
            sol  = solveLP(testCase.model);
            rxn1 = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            [p, m, f] = traceFluxPath(testCase.model, sol.x, rxn1, rxn1, 'verbose', false);
            testCase.verifyEqual(f, 1.0);
            testCase.verifyNumElements(p, 1);
            testCase.verifyEmpty(m);
        end

        function compareFluxesReturnsResult(testCase)
            testCase.assumeSolver('solveLP');
            solA = solveLP(testCase.model);
            o2exch = find(strcmp(testCase.model.rxnNames, 'O2 exchange'), 1);
            modelAna = setParam(testCase.model, 'eq', testCase.model.rxns(o2exch), 0);
            solB = solveLP(modelAna);
            result = compareFluxes(testCase.model, solA.x, solB.x, 'verbose', false);
            testCase.verifyTrue(isfield(result, 'turnedOn'));
            testCase.verifyTrue(isfield(result, 'turnedOff'));
            testCase.verifyTrue(isfield(result, 'flipped'));
            testCase.verifyTrue(isfield(result, 'changed'));
            testCase.verifyClass(result.changed.rxn, 'cell');
            testCase.verifyNotEmpty(result.changed.rxn);
            testCase.verifyEqual(numel(result.changed.flux1), numel(result.changed.rxn));
            testCase.verifyEqual(numel(result.changed.absDelta), numel(result.changed.rxn));
            % sorted descending
            testCase.verifyTrue(all(diff(result.changed.absDelta) <= 0));
        end

        function compareFluxesIdenticalIsEmpty(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            result = compareFluxes(testCase.model, sol.x, sol.x, 'verbose', false);
            testCase.verifyEmpty(result.changed.rxn);
            testCase.verifyEmpty(result.turnedOn);
            testCase.verifyEmpty(result.turnedOff);
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
            % Bare call exercises the default method (ACHR).
            evalc('sols = randomSampling(testCase.model, 5);');
            testCase.verifyEqual(size(sols, 1), numel(testCase.model.rxns));
        end

        function randomSamplingDeterministicWithSeed(testCase)
            % goodRxns reuse is specific to the random-objective method.
            evalc('[~, gR] = randomSampling(testCase.model, 0, ''method'', ''randomObjective'', ''runParallel'', false);');
            evalc('s1 = randomSampling(testCase.model, 5, ''method'', ''randomObjective'', ''seed'', 42, ''runParallel'', false, ''goodRxns'', gR);');
            evalc('s2 = randomSampling(testCase.model, 5, ''method'', ''randomObjective'', ''seed'', 42, ''runParallel'', false, ''goodRxns'', gR);');
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

        function getMinimalMediumReturnsMedium(testCase)
            testCase.assumeMILPSolver();
            evalc('[med, idx] = getMinimalMedium(testCase.model, ''verbose'', false);');
            testCase.verifyClass(med, 'cell');
            testCase.verifyClass(idx, 'double');
            % Every returned reaction must have lb < 0 (it is an uptake exchange)
            [~, exchIdx] = getExchangeRxns(testCase.model);
            testCase.verifyTrue(all(ismember(idx, exchIdx)));
            testCase.verifyTrue(all(testCase.model.lb(idx) < 0));
        end

        function getMinimalMediumExplicitGrowth(testCase)
            testCase.assumeMILPSolver();
            sol = solveLP(testCase.model);
            minG = 0.05 * sol.f;
            evalc('[med, idx] = getMinimalMedium(testCase.model, ''minGrowth'', minG, ''verbose'', false);');
            testCase.verifyClass(med, 'cell');
            testCase.verifyNotEmpty(med);
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

    methods (Static, Access = private)
        function m = loopModel()
            % R_in -> a; R1: a->b; R2: b->a (futile loop); R_out: b ->.
            m = struct();
            m.id='loop'; m.comps={'c'}; m.compNames={'cyt'};
            m.mets={'a';'b'}; m.metNames={'a';'b'}; m.metComps=[1;1];
            m.S = sparse([ 1 -1  1  0;
                           0  1 -1 -1]);
            m.rxns={'R_in';'R1';'R2';'R_out'}; m.rxnNames=m.rxns;
            m.lb=[0;0;0;0]; m.ub=[10;1000;1000;1000]; m.rev=[0;0;0;0];
            m.c=[0;0;0;1]; m.b=zeros(2,1);
            m.genes={}; m.grRules={'';'';'';''}; m.rxnGeneMat=sparse(4,0);
            m.metFormulas={'C';'C'};
        end
    end
end
