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

        function ravenListFormatsItems(testCase)
            msg = ravenList('Bad genes:', {'G1','G2'});
            testCase.verifySubstring(msg, 'G1');
            testCase.verifySubstring(msg, 'G2');
        end

        function ravenListTrimsToTenByDefault(testCase)
            items = arrayfun(@(n) sprintf('G%d',n), 1:12, 'UniformOutput', false);
            msg = ravenList('Too many:', items);
            testCase.verifySubstring(msg, '...and 3 more');
        end

        function ravenListNoTrimWhenFalse(testCase)
            items = arrayfun(@(n) sprintf('G%d',n), 1:12, 'UniformOutput', false);
            msg = ravenList('All:', items, false);
            testCase.verifySubstring(msg, 'G12');
            testCase.verifyFalse(contains(msg, '...and'));
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

        function parallelWorkersRAVENFalseReturnsZero(testCase)
            testCase.verifyEqual(parallelWorkersRAVEN(false), 0);
        end

        function parallelWorkersRAVENTrueNoPCTReturnsZero(testCase)
            testCase.assumeFalse(license('test','Distrib_Computing_Toolbox'), ...
                'Parallel Computing Toolbox present; skipping no-PCT path.');
            testCase.verifyEqual(parallelWorkersRAVEN(true), 0);
        end

        function parallelWorkersRAVENTruePCTReturnsInf(testCase)
            testCase.assumeTrue(license('test','Distrib_Computing_Toolbox'), ...
                'Parallel Computing Toolbox not present; skipping PCT path.');
            testCase.verifyEqual(parallelWorkersRAVEN(true), Inf);
        end

        function runRAVENtestsExists(testCase)
            % Do not execute it (would run the legacy suite); just confirm presence.
            testCase.verifyEqual(exist('runRAVENtests','file'), 2);
        end

        function parseGrRuleBuildsTree(testCase)
            t = parseGrRule('(G1 and G2) or G3');
            testCase.verifyEqual(t.type, 'or');
            testCase.verifyEqual(numel(t.children), 2);
            testCase.verifyEqual(t.children{1}.type, 'and');
            testCase.verifyEqual(t.children{2}.id, 'G3');
        end

        function parseGrRuleEmptyIsEmpty(testCase)
            testCase.verifyEmpty(parseGrRule(''));
            testCase.verifyEmpty(parseGrRule('   '));
        end

        function parseGrRuleTakesWholeWordOperators(testCase)
            % Gene IDs containing "and"/"or" are genes, not operators.
            t = parseGrRule('RAND1 or ORF2');
            testCase.verifyEqual(t.type, 'or');
            testCase.verifyEqual(t.children{1}.id, 'RAND1');
            testCase.verifyEqual(t.children{2}.id, 'ORF2');
        end

        function parseGrRuleAcceptsSymbolOperators(testCase)
            testCase.verifyEqual(parseGrRule('G1 & G2').type, 'and');
            testCase.verifyEqual(parseGrRule('G1 | G2').type, 'or');
            testCase.verifyEqual(parseGrRule('G1 AND G2').type, 'and');
        end

        function parseGrRuleRejectsMalformed(testCase)
            testCase.verifyError(@() parseGrRule('(G1 and G2'), 'RAVEN:badGrRule');
            testCase.verifyError(@() parseGrRule('G1 and'), 'RAVEN:badGrRule');
            testCase.verifyError(@() parseGrRule('G1 G2'), 'RAVEN:badGrRule');
        end

        function isDnfGrRuleClassifies(testCase)
            testCase.verifyTrue(isDnfGrRule('(G1 and G2) or G3'));
            testCase.verifyTrue(isDnfGrRule('G1'));
            testCase.verifyTrue(isDnfGrRule(''));
            testCase.verifyFalse(isDnfGrRule('(G1 or G2) and G3'));
            testCase.verifyFalse(isDnfGrRule('G1 and (G2 or G3)'));
        end

        function isDnfGrRuleIgnoresRedundantBrackets(testCase)
            % Bracketing alone never makes a rule non-DNF. Both of these were
            % false positives for the one-level split this replaced.
            testCase.verifyTrue(isDnfGrRule('((G1 and G2) or G3)'));
            testCase.verifyTrue(isDnfGrRule('(G1 or G2) or (G3 and G4)'));
        end

        function grRuleToDNFDistributes(testCase)
            testCase.verifyEqual(grRuleToDNF('g1 and (g2 or g3)'), ...
                {{'g1','g2'}, {'g1','g3'}});
            % Two or:s, but four isozymes: the count of ' or ' is not the
            % count of clauses.
            testCase.verifyEqual(numel(grRuleToDNF('(g1 or g2) and (g3 or g4)')), 4);
        end

        function grRuleToDNFSimpleCases(testCase)
            testCase.verifyEqual(grRuleToDNF(''), {});
            testCase.verifyEqual(grRuleToDNF('G1'), {{'G1'}});
            testCase.verifyEqual(grRuleToDNF('(G1 and G2) or G3'), {{'G1','G2'},{'G3'}});
        end

        function grRuleToStringRoundTrips(testCase)
            rules = {'(G1 and G2) or G3', 'G1', 'G1 and G2', 'G1 or G2'};
            for k = 1:numel(rules)
                testCase.verifyEqual(grRuleToString(parseGrRule(rules{k})), rules{k});
            end
            testCase.verifyEqual(grRuleToString(parseGrRule('')), '');
        end

    end
end
