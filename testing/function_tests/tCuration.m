classdef tCuration < RavenTestCase
% tCuration  Tests for the table-driven curation functions in curation/.

    methods (Test)

        function curateModelFromTablesNeedsTableFiles(testCase)
            % curateModelFromTables reads its mets/genes/rxns input from files
            % (tab-delimited / Excel); without those it cannot run.
            testCase.assumeFail('Requires curation input table files.');
        end

    end
end
