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
            loaded = load(fullfile(testCase.ravenRoot,'testing','unit_tests', ...
                'test_data','ecoli_textbook.mat'),'model');
            testCase.modelTemplate = loaded.model;
            try
                setRavenSolver('glpk');
            catch
                % leave whatever solver is configured; solver tests assume below
            end
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
    end
end
