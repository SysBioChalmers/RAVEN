classdef tIO < RavenTestCase
% tIO  Tests for the import/export and file-utility functions in io/.
%
%   Format support is exercised mainly via export->import round-trips on the
%   test model and via importing the tutorial 'empty' files.

    methods (TestClassSetup)
        function javaPaths(~)
            % Excel I/O relies on Apache POI being on the static Java path.
            try, addJavaPaths(); catch, end %#ok<NOCOM>
        end
    end

    methods (Test)

        function addJavaPathsRuns(testCase)
            testCase.verifyWarningFree(@() addJavaPaths());
        end

        function getFullPathReturnsAbsolute(testCase)
            p = getFullPath('.');
            if ispc
                % Windows: drive letter (C:\) or UNC (\\server).
                testCase.verifyTrue(~isempty(regexp(p, '^([A-Za-z]:|\\\\)', 'once')));
            else
                % Unix/macOS: absolute paths start with /.
                testCase.verifyTrue(startsWith(p, '/'));
            end
        end

        function getMD5HashReturnsHexDigest(testCase)
            f = fullfile(testCase.ravenRoot,'testing','unit_tests', ...
                'test_data','ecoli_textbook.mat');
            h = getMD5Hash(f);
            testCase.verifyTrue(ischar(h) || isstring(h));
            testCase.verifyMatches(char(h), '^[0-9a-fA-F]{32}$');
        end

        function getToolboxVersionRuns(testCase)
            v = getToolboxVersion('RAVEN', 'ravenCobraWrapper.m');
            testCase.verifyNotEmpty(v);
        end

        function importModelReadsSBML(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xml');
            evalc('m = importModel(f);');
            testCase.verifyTrue(isfield(m, 'rxns'));
        end

        function importExcelModelReadsXlsx(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xlsx');
            evalc('m = importExcelModel(f);');
            testCase.verifyTrue(isfield(m, 'rxns'));
        end

        function readYAMLmodelReadsYml(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.yml');
            evalc('m = readYAMLmodel(f);');
            testCase.verifyTrue(isfield(m, 'rxns'));
        end

        function parseYAMLReturnsTree(testCase)
            hasPyyaml = false;
            try, py.importlib.import_module('yaml'); hasPyyaml = true; catch, end %#ok<NOCOM>
            testCase.assumeDependency(hasPyyaml, 'Python pyyaml');
            f = fullfile(testCase.ravenRoot,'tutorial','empty.yml');
            out = parseYAML(f);
            testCase.verifyNotEmpty(out);
        end

        function exportImportSBMLRoundTrip(testCase)
            f = [tempname '.xml'];
            testCase.addTeardown(@() delete(f));
            evalc('exportModel(testCase.model, f);');
            evalc('m2 = importModel(f);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function exportImportExcelRoundTrip(testCase)
            f = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f));
            evalc('exportToExcelFormat(testCase.model, f);');
            evalc('m2 = importExcelModel(f);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function writeReadYAMLRoundTrip(testCase)
            f = [tempname '.yml'];
            testCase.addTeardown(@() delete(f));
            evalc('writeYAMLmodel(testCase.model, f);');
            evalc('m2 = readYAMLmodel(f);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function exportToTabDelimitedWritesFiles(testCase)
            d = [tempname filesep];
            mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            evalc('exportToTabDelimited(testCase.model, d);');
            testCase.verifyTrue(exist(fullfile(d,'excelRxns.txt'),'file')==2);
        end

        function exportModelToSIFWritesFile(testCase)
            f = [tempname '.sif'];
            testCase.addTeardown(@() delete(f));
            evalc('exportModelToSIF(testCase.model, f);');
            testCase.verifyTrue(exist(f,'file')==2);
        end

        function exportForGitWritesRepo(testCase)
            d = [tempname filesep];
            mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            evalc('exportForGit(testCase.model, ''ec'', d, {''yml''});');
            testCase.verifyNotEmpty(dir(fullfile(d,'**','*.yml')));
        end

        function loadWorkbookLoadsFile(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xlsx');
            wb = loadWorkbook(f);
            testCase.verifyNotEmpty(wb);
        end

        function checkFileExistenceFindsFile(testCase)
            f = fullfile(testCase.ravenRoot,'testing','unit_tests', ...
                'test_data','ecoli_textbook.mat');
            out = checkFileExistence(f, 1, false, true);
            testCase.verifyNotEmpty(out);
        end

        function cleanSheetTrimsComments(testCase)
            raw = {'#comment', '#comment'; 'a', 'b'; '', ''};
            out = cleanSheet(raw);
            testCase.verifyClass(out, 'cell');
        end

        function loadSheetReadsSheet(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xlsx');
            wb = loadWorkbook(f);
            [raw, flag] = loadSheet(wb, 'RXNS');
            testCase.verifyEqual(flag, 0);
            testCase.verifyClass(raw, 'cell');
        end

        function writeSheetWritesToWorkbook(testCase)
            f = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f));
            wb = loadWorkbook(f, true);   % create empty
            % writeSheet only knows column widths for the standard RAVEN sheets.
            wb = writeSheet(wb, 'RXNS', 0, {'header'}, [], {'value'});
            testCase.verifyNotEmpty(wb);
        end

    end
end
