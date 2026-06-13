classdef (Abstract) RavenTestCase < matlab.unittest.TestCase
% RavenTestCase  Shared base class for the RAVEN function test suite.
%
%   Provides, once per test class:
%     - ravenRoot : the RAVEN installation root.
%     - model     : a small E. coli textbook model in RAVEN format, reloaded
%                   fresh for every test method so mutations cannot leak
%                   between tests.
%     - a deterministic LP solver (GLPK, which ships with RAVEN).
%
%   It also offers assume* helpers that mark a test as filtered (not failed)
%   when an optional dependency is missing, so the suite stays green in any
%   environment. Concrete test classes subclass this and add a methods(Test)
%   block.

    properties
        ravenRoot char
        model struct          % reset before every test method
    end

    properties (Access = protected)
        modelTemplate struct  % pristine copy, loaded once per class
    end

    methods (TestClassSetup)
        function setupRaven(testCase)
            testCase.ravenRoot = findRAVENroot();
            loaded = load(fullfile(testCase.ravenRoot,'testing','function_tests', ...
                'test_data','ecoli_textbook.mat'),'model');
            testCase.modelTemplate = loaded.model;
            try
                setRavenSolver('glpk');
            catch
                % leave whatever solver is configured; solver tests assume below
            end

            % Silence progress reporting for the duration of the test class,
            % restoring whatever the user had configured afterwards.
            origProgress = getpref('RAVEN','progressBar','auto');
            setpref('RAVEN','progressBar','none');
            testCase.addTeardown(@() setpref('RAVEN','progressBar',origProgress));
        end

        function suppressFigures(testCase)
            % Keep plotting functions from popping windows or leaking figures.
            orig = get(groot, 'defaultFigureVisible');
            set(groot, 'defaultFigureVisible', 'off');
            testCase.addTeardown(@() set(groot, 'defaultFigureVisible', orig));
            testCase.addTeardown(@() close('all'));
        end
    end

    methods (TestMethodSetup)
        function freshModel(testCase)
            % Hand each test an untouched copy of the model.
            testCase.model = testCase.modelTemplate;
        end
    end

    methods (Access = protected)
        function assumeSolver(testCase, solverName)
            % Skip unless the named solver is on the MATLAB path.
            testCase.assumeNotEmpty(which(solverName), ...
                sprintf('Solver "%s" not available; test skipped.', solverName));
        end

        function assumeDependency(testCase, isAvailable, label)
            % Skip when an external binary / network / GUI dependency is absent.
            testCase.assumeTrue(logical(isAvailable), ...
                sprintf('Dependency "%s" not available; test skipped.', label));
        end

        function assumeMILPSolver(testCase)
            % Switch to a MILP-capable solver for this test, or skip if none.
            % GLPK (the default) cannot solve MILPs in RAVEN.
            if ~isempty(which('gurobi'))
                evalc('setRavenSolver(''gurobi'')');
            else
                ok = false;
                try, evalc('setRavenSolver(''scip'')'); ok = true; catch, end %#ok<NOCOM>
                testCase.assumeTrue(ok, 'No MILP solver (gurobi/scip) available; test skipped.');
            end
            testCase.addTeardown(@() evalc('setRavenSolver(''glpk'')'));
        end

        function m = taskTestModel(~)
            % Small synthetic model used by task/gap-filling tests.
            m = struct();
            m.id = 'testModel'; m.rxns = {}; m.S = []; m.rev = [];
            m.mets = {'as';'ac';'bc';'cc';'dc';'ec';'es';'fc'};
            m.metNames = {'a';'a';'b';'c';'d';'e';'e';'f'};
            m.comps = {'s';'c'}; m.compNames = m.comps;
            m.metComps = [1;2;2;2;2;2;1;2];
            m.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
            m.grRules = {}; m.rxnGeneMat = [];
            r = struct();
            r.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
            r.equations = {'=> a[s]';'a[s] <=> a[c]';'a[c] <=> b[c] + c[c]'; ...
                'a[c] <=> 2 d[c]';'b[c] + c[c] => e[c]';'2 d[c] => e[c]'; ...
                'e[c] => e[s]';'e[s] =>';'a[c] <=> f[c]';'f[c] <=> e[c]'};
            r.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};
            evalc('m = addRxns(m, r, 3);');
            m.c = [0;0;0;0;0;0;0;1;0;0];
            m.ub = repmat(1000,10,1);
            m.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
            m.rxnNames = m.rxns;
            m.b = repmat(0,8,1);
            m = closeModel(m);
        end

        function t = taskTestStruct(~)
            % Task that the model above can satisfy: make e[s] from a[s].
            t = struct();
            t.id = 'Gen e[s] from a[s]'; t.description = t.id;
            t.shouldFail = false; t.printFluxes = false; t.comments = '';
            t.inputs = {'a[s]'}; t.LBin = 0; t.UBin = inf;
            t.outputs = {'e[s]'}; t.LBout = 1; t.UBout = 1;
            t.equations = {}; t.LBequ = []; t.UBequ = [];
            t.changed = {}; t.LBrxn = {}; t.UBrxn = {};
        end
    end
end
