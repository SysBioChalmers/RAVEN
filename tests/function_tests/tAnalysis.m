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

        function haveFluxSeedIsReproducible(testCase)
            % The order reactions are tested in is randomised; a given
            % seed must make that order (and so results) reproducible.
            I1 = haveFlux(testCase.model, 'seed', 42);
            I2 = haveFlux(testCase.model, 'seed', 42);
            testCase.verifyEqual(I1, I2);
        end

        function getMinNrFluxesReturnsFlux(testCase)
            testCase.assumeMILPSolver();
            evalc('[x, I, exitFlag] = getMinNrFluxes(testCase.model, testCase.model.rxns);');
            testCase.verifyNotEmpty(x);
            %A problem this small is solved well within the time limit
            testCase.verifyEqual(exitFlag, 1);
        end

        function getMinNrFluxesFormulationsAgree(testCase)
            % The two formulations build different MILPs for the same
            % problem, so they may pick different reactions, but both
            % minimise the same objective and must therefore pick the same
            % number of them.
            testCase.assumeMILPSolver();
            evalc('[xIrrev, Iirrev, flagIrrev] = getMinNrFluxes(testCase.model, testCase.model.rxns, [], [], ''irrev'');');
            evalc('[xRev, Irev, flagRev] = getMinNrFluxes(testCase.model, testCase.model.rxns, [], [], ''reversible'');');
            testCase.verifyEqual(flagIrrev, 1);
            testCase.verifyEqual(flagRev, 1);
            testCase.verifyNumElements(xRev, numel(testCase.model.rxns));
            testCase.verifyNumElements(Irev, numel(testCase.model.rxns));
            testCase.verifyEqual(sum(Irev), sum(Iirrev));
        end

        function getMinNrFluxesResolveTiesIsDeterministic(testCase)
            % The reversible-formulation MILP is degenerate on a real model
            % (see getMinNrFluxesFormulationsAgree above): resolveTies pins
            % it to a canonical answer instead of leaving the choice to the
            % solver/seed, so three independent runs must all agree
            % (raven-gecko-parity#104).
            testCase.assumeMILPSolver();
            evalc(['[~, I1, flag1] = getMinNrFluxes(testCase.model, testCase.model.rxns, ' ...
                '[], [], ''reversible'', false, ''resolveTies'', true);']);
            evalc(['[~, I2, flag2] = getMinNrFluxes(testCase.model, testCase.model.rxns, ' ...
                '[], [], ''reversible'', false, ''resolveTies'', true);']);
            evalc(['[~, I3, flag3] = getMinNrFluxes(testCase.model, testCase.model.rxns, ' ...
                '[], [], ''reversible'', false, ''resolveTies'', true);']);
            testCase.verifyEqual(flag1, 1);
            testCase.verifyEqual(flag2, 1);
            testCase.verifyEqual(flag3, 1);
            testCase.verifyEqual(I1, I2);
            testCase.verifyEqual(I2, I3);
        end

        function getMinNrFluxesResolveTiesOnlySupportsReversible(testCase)
            testCase.verifyError( ...
                @() getMinNrFluxes(testCase.model, testCase.model.rxns, [], [], 'irrev', ...
                    false, 'resolveTies', true), ...
                'RAVEN:badInput');
        end

        function getMinNrFluxesRejectsUnknownFormulation(testCase)
            testCase.verifyError( ...
                @() getMinNrFluxes(testCase.model, testCase.model.rxns, [], [], 'both'), ...
                'RAVEN:badInput');
        end

        function getAllSubGraphsReturnsResult(testCase)
            sg = getAllSubGraphs(testCase.model);
            testCase.verifyNotEmpty(sg);
        end

        function findGeneDeletionsRuns(testCase)
            evalc('[genes, fluxes] = findGeneDeletions(testCase.model, ''sgd'');');
            testCase.verifyNotEmpty(genes);
        end

        function findGeneDeletionsComputesGrRatio(testCase)
            evalc(['[genes, ~, ~, ~, grRatioMuts] = findGeneDeletions(testCase.model, ' ...
                '''sgd'');']);
            testCase.verifyNotEmpty(genes);
            testCase.verifyEqual(numel(grRatioMuts), numel(genes));
            testCase.verifyGreaterThan(max(grRatioMuts), 0);
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

        function walkFluxesRefusesNonInteractive(testCase)
            % The test suite itself runs under -batch, so this exercises the
            % real guard rather than a mocked one.
            testCase.assumeTrue(batchStartupOptionUsed, ...
                'Test runner is not using -batch; the non-interactive guard cannot be exercised here.');
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            testCase.verifyError(@() walkFluxes(testCase.model, sol.x, biomassRxn), ...
                'walkFluxes:nonInteractive');
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

        function compareFluxesMetaboliteListRestrictsResult(testCase)
            testCase.assumeSolver('solveLP');
            solA = solveLP(testCase.model);
            o2exch = find(strcmp(testCase.model.rxnNames, 'O2 exchange'), 1);
            modelAna = setParam(testCase.model, 'eq', testCase.model.rxns(o2exch), 0);
            solB = solveLP(modelAna);
            unfiltered = compareFluxes(testCase.model, solA.x, solB.x, 'verbose', false);
            testCase.assumeNotEmpty(unfiltered.changed.rxn);

            % Name a metabolite of the largest-changing reaction, so the
            % filtered result is guaranteed to be non-empty and the subset
            % checks below are not satisfied vacuously.
            topRxn = strcmp(testCase.model.rxns, unfiltered.changed.rxn{1});
            met = testCase.model.metNames{find(testCase.model.S(:,topRxn) ~= 0, 1)};
            filtered = compareFluxes(testCase.model, solA.x, solB.x, ...
                'metaboliteList', {met}, 'verbose', false);
            testCase.verifyNotEmpty(filtered.changed.rxn);

            % Every kept reaction must involve the named metabolite, and the
            % filtered result can only be a subset of the unfiltered one.
            metIdx = strcmpi(met, testCase.model.metNames);
            withMet = testCase.model.rxns(any(testCase.model.S(metIdx,:) ~= 0, 1));
            testCase.verifyTrue(all(ismember(filtered.changed.rxn, withMet)));
            testCase.verifyTrue(all(ismember(filtered.changed.rxn, unfiltered.changed.rxn)));
        end

        function compareFluxesUnknownMetaboliteWarns(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            testCase.verifyWarning(@() compareFluxes(testCase.model, sol.x, sol.x, ...
                'metaboliteList', {'no such metabolite'}, 'verbose', false), ...
                'RAVEN:warning');
        end

        function getFluxZZeroVarianceSignMatchesGeneralCase(testCase)
            % The zero-variance branch's sign must agree with the general
            % branch: positive when flux increased from A to B, negative
            % when it decreased.
            solA = [1 1 1; 5 5 5];    % rxn1 constant at 1, rxn2 constant at 5
            solB = [5 5 5; 1 1 1];    % rxn1 increased to 5, rxn2 decreased to 1
            Z = getFluxZ(solA, solB);
            testCase.verifyEqual(Z(1), 100);
            testCase.verifyEqual(Z(2), -100);
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

        function reporterMetabolitesKeepsPercentInMetNames(testCase)
            % A metNames entry containing "%" must survive intact in the
            % outputFile report, not be truncated by fprintf misreading it
            % as a format directive.
            m = testCase.model;
            m.metNames{1} = 'metabolite 30% pure';
            pvals = rand(numel(m.genes), 1);
            outFile = [tempname '.txt'];
            c = onCleanup(@() delete(outFile));
            evalc('reporterMetabolites(m, m.genes, pvals, ''outputFile'', outFile);');
            content = fileread(outFile);
            testCase.verifySubstring(content, 'metabolite 30% pure');
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
end
