classdef tSampling < RavenTestCase
% tSampling  Tests for the MCMC flux samplers (CHRR, ACHR, MVE) dispatched
% through randomSampling.
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

        % ---- sampler implementations return diagnostics -----------------

        function sampleCHRRReturnsInfo(testCase)
            evalc('[sols, sInfo] = sampleCHRR(testCase.model, 10, 5, 20, 1);');
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            testCase.verifyTrue(isfield(sInfo, 'nDimensions'));
            testCase.verifyTrue(isfield(sInfo, 'mveConverged'));
            % Samples must satisfy steady state and bounds.
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
            testCase.verifyGreaterThanOrEqual(min(sols, [], 2), testCase.model.lb - 1e-6);
            testCase.verifyLessThanOrEqual(max(sols, [], 2), testCase.model.ub + 1e-6);
        end

        function sampleACHRReturnsSteadyState(testCase)
            evalc('sols = sampleACHR(testCase.model, 10, 5, [], 1);');
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
        end

        % ---- randomSampling dispatch ------------------------------------

        function randomSamplingDefaultsToACHR(testCase)
            % The default method must be ACHR: with a fixed seed the default call
            % must equal an explicit 'achr' call.
            evalc('sDefault = randomSampling(testCase.model, 8, ''thinning'', 5, ''seed'', 3);');
            evalc('sAchr = randomSampling(testCase.model, 8, ''method'', ''achr'', ''thinning'', 5, ''seed'', 3);');
            testCase.verifyEqual(sDefault, sAchr);
        end

        function randomSamplingCHRRDispatch(testCase)
            evalc(['sols = randomSampling(testCase.model, 10, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
        end

        function randomSamplingCHRRExposesInfo(testCase)
            % The CHRR convergence diagnostics must be reachable through the
            % documented randomSampling entry point, not only sampleCHRR.
            evalc(['[sols, gr, sInfo] = randomSampling(testCase.model, 10, ' ...
                '''method'', ''chrr'', ''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            testCase.verifyTrue(isfield(sInfo, 'mveConverged'));
        end

        function randomSamplingCHRRDeterministic(testCase)
            evalc(['s1 = randomSampling(testCase.model, 8, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            evalc(['s2 = randomSampling(testCase.model, 8, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            testCase.verifyEqual(s1, s2);
        end

        function randomSamplingCHRRVaries(testCase)
            % A correct CHRR run must mix: at least one reaction varies, and two
            % seeds give distinct sample matrices. Guards against a degenerate
            % sampler that collapses to a single point (which would still satisfy
            % S*v=0 and the bounds, and remain deterministic).
            evalc(['s1 = randomSampling(testCase.model, 50, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            evalc(['s2 = randomSampling(testCase.model, 50, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 2);']);
            testCase.verifyGreaterThan(max(std(s1, 0, 2)), 1e-3);
            testCase.verifyNotEqual(s1, s2);
        end

        function randomSamplingUnknownMethodErrors(testCase)
            testCase.verifyError( ...
                @() randomSampling(testCase.model, 5, 'method', 'gibbs'), ...
                'RAVEN:badInput');
        end

        function randomSamplingLoopDetectionUsesModelBounds(testCase)
            % Loop detection must take its threshold from the model. R2/R3
            % form an a -> b -> a loop that runs up to this model's 100 cap,
            % while the linear path R1 -> R4 is held to 10. A threshold that
            % assumes ±1000 bounds matches neither, leaving the loop
            % reactions to be sampled as if they were loop-free.
            m = tSampling.loopModel();
            evalc(['[~, goodRxns] = randomSampling(m, 2, ''method'', ' ...
                '''randomObjective'', ''seed'', 1);']);
            testCase.verifyEqual(sort(goodRxns(:))', [1 4]);
        end

    end

    methods (Static, Access = private)

        function m = loopModel()
            m = struct();
            m.id       = 'loop';
            m.rxns     = {'R1';'R2';'R3';'R4'};
            m.rxnNames = m.rxns;
            m.mets     = {'a';'b'};
            m.metNames = {'a';'b'};
            m.metComps = [1;1];
            m.comps    = {'c'};
            m.compNames= {'cytosol'};
            %             R1  R2  R3  R4
            m.S = sparse([ 1  -1   1   0;    % a
                           0   1  -1  -1]);  % b
            m.lb  = [0;0;0;0];
            m.ub  = [10;100;100;100];
            m.rev = [0;0;0;0];
            m.c   = [0;0;0;1];
            m.b   = zeros(2,1);
            m.genes = {}; m.grRules = {'';'';'';''}; m.rxnGeneMat = sparse(4,0);
            m.metFormulas = {'C';'C'};
        end

    end
end
