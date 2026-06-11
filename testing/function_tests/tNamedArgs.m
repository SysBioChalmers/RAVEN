classdef tNamedArgs < RavenTestCase
% tNamedArgs  Verifies that functions converted to parseRAVENargs accept both
%             the positional and the name-value calling conventions, and that
%             the two forms give identical results.

    methods (Test)

        function getIndexesPositionalAndNamedMatch(testCase)
            pos = getIndexes(testCase.model, testCase.model.rxns(1:3), 'rxns', true);
            nv  = getIndexes(testCase.model, testCase.model.rxns(1:3), 'rxns', ...
                'returnLogical', true);
            testCase.verifyEqual(nv, pos);
        end

        function addTransportPositionalAndNamedMatch(testCase)
            evalc(['m1 = addTransport(testCase.model, ''c'', ''e'', ' ...
                '{''6-phospho-D-glucono-1,5-lactone''}, false, false, ''tr_'');']);
            evalc(['m2 = addTransport(testCase.model, ''c'', ''e'', ' ...
                '''metNames'', {''6-phospho-D-glucono-1,5-lactone''}, ' ...
                '''isRev'', false, ''onlyToExisting'', false, ''prefix'', ''tr_'');']);
            testCase.verifyEqual(numel(m2.rxns), numel(m1.rxns));
        end

        function addTransportNamedSubsetOfOptions(testCase)
            % A subset of options by name; the rest keep their defaults.
            evalc(['m = addTransport(testCase.model, ''c'', ''e'', ' ...
                '''metNames'', {''6-phospho-D-glucono-1,5-lactone''}, ' ...
                '''onlyToExisting'', false);']);
            testCase.verifyGreaterThan(numel(m.rxns), numel(testCase.model.rxns));
        end

        function exportToExcelFormatPositionalAndNamed(testCase)
            f1 = [tempname '.xlsx'];
            f2 = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f1));
            testCase.addTeardown(@() delete(f2));
            try, addJavaPaths(); catch, end %#ok<NOCOM>
            evalc('exportToExcelFormat(testCase.model, f1, true);');
            evalc('exportToExcelFormat(testCase.model, ''fileName'', f2, ''sortIds'', true);');
            testCase.verifyTrue(exist(f1,'file')==2);
            testCase.verifyTrue(exist(f2,'file')==2);
        end

        function tooManyPositionalErrors(testCase)
            testCase.verifyError(@() getIndexes(testCase.model, ...
                testCase.model.rxns(1), 'rxns', true, 'extra'), ?MException);
        end

    end
end
