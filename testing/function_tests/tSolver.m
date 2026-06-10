classdef tSolver < RavenTestCase
% tSolver  Tests for the LP/QP solving and solver-configuration functions in solver/.

    methods (Test)

        function setRavenSolverAcceptsGlpk(testCase)
            testCase.verifyWarningFree(@() setRavenSolver('glpk'));
        end

        function setRavenSolverRejectsUnknown(testCase)
            testCase.verifyError(@() setRavenSolver('notarealsolver'), ?MException);
        end

        function solveLPFindsOptimum(testCase)
            sol = solveLP(testCase.model);
            testCase.verifyEqual(sol.stat, 1);
            testCase.verifyNumElements(sol.x, numel(testCase.model.rxns));
            testCase.verifyEqual(sol.f, 0.87392, 'AbsTol', 1e-3);
        end

        function solveLPMinimiseFluxes(testCase)
            sol = solveLP(testCase.model, 1);
            testCase.verifyEqual(sol.stat, 1);
        end

        function optimizeProbSolvesLP(testCase)
            res = optimizeProb(testCase.modelLP());
            testCase.verifyClass(res, 'struct');
            testCase.verifyTrue(isfield(res, 'obj'));
        end

        function checkSolutionOnFeasibleProblem(testCase)
            res = optimizeProb(testCase.modelLP());
            [isFeasible, isOptimal] = checkSolution(res);
            testCase.verifyTrue(logical(isFeasible));
            testCase.verifyTrue(logical(isOptimal));
        end

        function splitProbForConditioningReturnsProb(testCase)
            [p2, condInfo] = splitProbForConditioning(testCase.modelLP(), 1e6); %#ok<ASGLU>
            testCase.verifyClass(p2, 'struct');
        end

        function qMOMAReturnsFluxVectors(testCase)
            testCase.assumeDependency(exist('quadprog','file')==2, ...
                'Optimization Toolbox (quadprog)');
            model2 = setParam(testCase.model, 'eq', testCase.model.rxns(1), 0);
            evalc('[fA, fB, flag] = qMOMA(testCase.model, model2);');
            testCase.verifyNumElements(fA, numel(testCase.model.rxns));
            testCase.verifyNumElements(fB, numel(testCase.model.rxns));
        end

    end

    methods (Access = private)
        function prob = modelLP(testCase)
            % Build the COBRA-style LP problem for the test model exactly as
            % solveLP does, so the low-level solver functions get a valid input.
            m = testCase.model;
            nMets = size(m.S, 1);
            prob = [];
            prob.c = [m.c * -1; zeros(nMets, 1)];
            prob.b = zeros(nMets, 1);
            prob.lb = [m.lb; m.b(:, 1)];
            prob.ub = [m.ub; m.b(:, min(size(m.b, 2), 2))];
            prob.csense = repmat('E', 1, numel(prob.b));
            prob.osense = 1;
            prob.vartype = repmat('C', 1, numel(prob.c));
            prob.A = [m.S, -speye(nMets)];
            prob.a = m.S;
        end
    end
end
