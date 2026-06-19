classdef tSampling < RavenTestCase
% tSampling  Tests for the MCMC flux samplers in analysis/ (CHRR, ACHR, MVE).
%
%   The maximum-volume ellipsoid solver is the numerical crux of CHRR, so it is
%   validated directly against analytic cases (box -> unit ball, scaled/sheared
%   box, triangle -> Steiner inellipse) before the end-to-end samplers are
%   exercised on the textbook model. All methods are LP-only, so the default
%   GLPK solver suffices (no MILP assumption needed).

    methods (Test)

        % ---- sampleMaxVolEllipse: analytic validation -------------------

        function sampleMaxVolEllipseUnitBox(testCase)
            % Box {-1 <= x_i <= 1} -> unit ball: centre 0, E*E' = I.
            A = [eye(3); -eye(3)];
            b = ones(6, 1);
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(3, 1));
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - eye(3)), 1e-3);
        end

        function sampleMaxVolEllipseScaledBox(testCase)
            % Box {-r_i <= x_i <= r_i} -> E*E' = diag(r_i^2).
            r = [2; 0.5; 3];
            A = [eye(3); -eye(3)];
            b = [r; r];
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(3, 1));
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - diag(r.^2)), 1e-3);
        end

        function sampleMaxVolEllipseShearedBox(testCase)
            % {-1 <= M*x <= 1} -> E*E' = inv(M)*inv(M)'.
            M = [1 0.5; 0 1];
            A = [M; -M];
            b = ones(4, 1);
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(2, 1));
            Minv = inv(M);
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - Minv*Minv'), 1e-3);
        end

        function sampleMaxVolEllipseTriangle(testCase)
            % Standard simplex {x>=0, x1+x2<=1}: the max-area inscribed ellipse is
            % the Steiner inellipse, centred at the centroid (1/3,1/3) with
            % E2 = [[1/9,-1/18],[-1/18,1/9]] (det 1/108) -- NOT the incircle. This
            % confirms the solver maximises volume, not just the largest ball.
            A = [-1 0; 0 -1; 1 1];
            b = [0; 0; 1];
            [x, E, converged] = sampleMaxVolEllipse(A, b, [0.25; 0.25]);
            E2exp = [1/9 -1/18; -1/18 1/9];
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x - [1/3; 1/3]), 1e-3);
            testCase.verifyLessThan(norm(E*E' - E2exp), 1e-3);
            testCase.verifyLessThan(abs(det(E*E') - 1/108), 1e-4);
        end

        function sampleMaxVolEllipseRejectsExterior(testCase)
            A = [eye(2); -eye(2)];
            b = ones(4, 1);
            testCase.verifyError(@() sampleMaxVolEllipse(A, b, [5; 0]), ...
                'RAVEN:sampling');
        end

        % ---- sampleFluxSpace: CHRR --------------------------------------

        function sampleFluxSpaceCHRRRuns(testCase)
            % Note: name the info output 'sInfo' (not 'info') — 'info' is a
            % MATLAB built-in, and evalc-created variables are not statically
            % visible to the parser, so 'info' would resolve to the builtin.
            evalc(['[sols, sInfo] = sampleFluxSpace(testCase.model, ''method'', ''chrr'', ' ...
                '''nSamples'', 10, ''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            testCase.verifyClass(sInfo, 'struct');
            % Samples must satisfy steady state S*v = 0.
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
            % Samples must respect the bounds.
            testCase.verifyGreaterThanOrEqual(min(sols, [], 2), testCase.model.lb - 1e-6);
            testCase.verifyLessThanOrEqual(max(sols, [], 2), testCase.model.ub + 1e-6);
        end

        function sampleFluxSpaceCHRRDeterministic(testCase)
            evalc(['s1 = sampleFluxSpace(testCase.model, ''method'', ''chrr'', ' ...
                '''nSamples'', 8, ''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            evalc(['s2 = sampleFluxSpace(testCase.model, ''method'', ''chrr'', ' ...
                '''nSamples'', 8, ''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            testCase.verifyEqual(s1, s2);
        end

        function sampleFluxSpaceCHRRVaries(testCase)
            % A correct CHRR run must mix: at least one reaction's flux must vary
            % across samples, and two distinct seeds must give distinct sample
            % matrices. Guards against a degenerate sampler that collapses to a
            % single point — which would still satisfy S*v=0 (v0 alone is a steady
            % state) and the bounds, and remain deterministic, thereby passing
            % every other CHRR test.
            evalc(['s1 = sampleFluxSpace(testCase.model, ''method'', ''chrr'', ' ...
                '''nSamples'', 50, ''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            evalc(['s2 = sampleFluxSpace(testCase.model, ''method'', ''chrr'', ' ...
                '''nSamples'', 50, ''thinning'', 5, ''nBurnin'', 20, ''seed'', 2);']);
            testCase.verifyGreaterThan(max(std(s1, 0, 2)), 1e-3);
            testCase.verifyNotEqual(s1, s2);
        end

        % ---- sampleFluxSpace: ACHR --------------------------------------

        function sampleFluxSpaceACHRRuns(testCase)
            evalc(['[sols, sInfo] = sampleFluxSpace(testCase.model, ''method'', ''achr'', ' ...
                '''nSamples'', 10, ''thinning'', 5, ''seed'', 1);']);
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            testCase.verifyTrue(isfield(sInfo, 'nWarmup'));
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
        end

        function sampleFluxSpaceACHRDeterministic(testCase)
            evalc(['s1 = sampleFluxSpace(testCase.model, ''method'', ''achr'', ' ...
                '''nSamples'', 8, ''thinning'', 5, ''seed'', 7);']);
            evalc(['s2 = sampleFluxSpace(testCase.model, ''method'', ''achr'', ' ...
                '''nSamples'', 8, ''thinning'', 5, ''seed'', 7);']);
            testCase.verifyEqual(s1, s2);
        end

        % ---- dispatcher errors ------------------------------------------

        function sampleFluxSpaceUnknownMethodErrors(testCase)
            testCase.verifyError( ...
                @() sampleFluxSpace(testCase.model, 'method', 'gibbs'), ...
                'RAVEN:badInput');
        end

    end
end
