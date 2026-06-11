classdef tUtils < RavenTestCase
% tUtils  Tests for the low-level helpers in utils/.

    methods (Test)

        function convertCharArrayFromChar(testCase)
            testCase.verifyEqual(convertCharArray('abc'), {'abc'});
        end

        function convertCharArrayFromCell(testCase)
            testCase.verifyEqual(convertCharArray({'a','b'}), {'a','b'});
        end

        function convertCharArrayFromString(testCase)
            testCase.verifyEqual(convertCharArray("str"), {'str'});
        end

        function dispEMThrowsByDefault(testCase)
            testCase.verifyError(@() dispEM('boom'), ?MException);
        end

        function dispEMThrowsWhenRequested(testCase)
            testCase.verifyError(@() dispEM('boom', true), ?MException);
        end

        function dispEMWarningIsPrintedNotThrown(testCase)
            out = evalc('dispEM(''heads up'', false)');
            testCase.verifySubstring(out, 'WARNING');
        end

        function emptyOrLogicalScalarAcceptsValid(testCase)
            testCase.verifyWarningFree(@() emptyOrLogicalScalar(true));
            testCase.verifyWarningFree(@() emptyOrLogicalScalar([]));
        end

        function emptyOrLogicalScalarRejectsInvalid(testCase)
            testCase.verifyError(@() emptyOrLogicalScalar([true false]), ?MException);
        end

        function emptyOrTextScalarAcceptsValid(testCase)
            testCase.verifyWarningFree(@() emptyOrTextScalar('abc'));
            testCase.verifyWarningFree(@() emptyOrTextScalar([]));
        end

        function emptyOrTextScalarRejectsInvalid(testCase)
            testCase.verifyError(@() emptyOrTextScalar(5), ?MException);
        end

        function emptyOrTextOrCellOfTextAcceptsValid(testCase)
            testCase.verifyWarningFree(@() emptyOrTextOrCellOfText({'a','b'}));
            testCase.verifyWarningFree(@() emptyOrTextOrCellOfText('a'));
        end

        function emptyOrTextOrCellOfTextRejectsInvalid(testCase)
            testCase.verifyError(@() emptyOrTextOrCellOfText(5), ?MException);
        end

        function printOrangeReturnsTextContainingInput(testCase)
            evalc('s = printOrange(''hello'');');
            testCase.verifySubstring(s, 'hello');
        end

        function parallelPoolRAVENnoParallelRuns(testCase)
            testCase.assumeDependency(license('test','Distrib_Computing_Toolbox'), ...
                'Parallel Computing Toolbox');
            testCase.verifyWarningFree(@() parallelPoolRAVEN(false));
        end

        function runRAVENtestsExists(testCase)
            % Do not execute it (would run the legacy suite); just confirm presence.
            testCase.verifyEqual(exist('runRAVENtests','file'), 2);
        end

    end
end
