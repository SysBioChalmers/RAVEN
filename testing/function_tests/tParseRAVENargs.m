classdef tParseRAVENargs < matlab.unittest.TestCase
% tParseRAVENargs  Tests for the parseRAVENargs optional-argument resolver.

    methods (Test)

        function defaultsWhenEmpty(testCase)
            p = parseRAVENargs({}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 1);
            testCase.verifyEqual(p.b, 'x');
        end

        function positionalAssignment(testCase)
            p = parseRAVENargs({10, 'y'}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 10);
            testCase.verifyEqual(p.b, 'y');
        end

        function partialPositionalKeepsLaterDefaults(testCase)
            p = parseRAVENargs({10}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 10);
            testCase.verifyEqual(p.b, 'x');
        end

        function nameValueAssignment(testCase)
            p = parseRAVENargs({'b', 'y'}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 1);     % untouched -> default
            testCase.verifyEqual(p.b, 'y');
        end

        function nameValueIsCaseInsensitive(testCase)
            p = parseRAVENargs({'B', 'y'}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.b, 'y');
        end

        function nameValueBothParams(testCase)
            p = parseRAVENargs({'b', 'y', 'a', 5}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 5);
            testCase.verifyEqual(p.b, 'y');
        end

        function oddCountStaysPositional(testCase)
            % A single string value that happens to match a name is positional.
            p = parseRAVENargs({'a'}, {'a', 1; 'b', 'x'});
            testCase.verifyEqual(p.a, 'a');
        end

        function tooManyPositionalErrors(testCase)
            testCase.verifyError(@() parseRAVENargs({1, 2, 3}, {'a', 1; 'b', 2}), ...
                'parseRAVENargs:tooManyArgs');
        end

        function validatorRejectsInvalid(testCase)
            spec = {'flag', false, @emptyOrLogicalScalar};
            testCase.verifyError(@() parseRAVENargs({'flag', [1 2 3]}, spec), ?MException);
        end

        function validatorAcceptsValid(testCase)
            spec = {'flag', false, @emptyOrLogicalScalar};
            p = parseRAVENargs({'flag', true}, spec);
            testCase.verifyTrue(p.flag);
        end

    end
end
