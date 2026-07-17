classdef tIO < RavenTestCase
% tIO  Tests for the import/export and file-utility functions in io/.
%
%   Format support is exercised mainly via export->import round-trips on the
%   test model and via importing the tutorial 'empty' files.

    methods (Test)

        function exportForGitWritesDependencies(testCase)
            % The toolbox-version lookup is now embedded in exportForGit and
            % is exercised via the dependencies.txt it writes
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xml');
            evalc('model = importModel(f);');
            outDir = tempname; mkdir(outDir);
            c = onCleanup(@() rmdir(outDir,'s'));
            evalc("exportForGit(model,'path',outDir,'formats',{'xml'},'subDirs',false)");
            dep = fileread(fullfile(outDir,'dependencies.txt'));
            testCase.verifySubstring(dep, 'RAVEN_toolbox');
        end

        function importModelReadsSBML(testCase)
            f = fullfile(testCase.ravenRoot,'tutorial','empty.xml');
            evalc('m = importModel(f);');
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

        function exportToExcelFormatWritesFile(testCase)
            f = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f));
            evalc('exportToExcelFormat(testCase.model, f);');
            testCase.verifyTrue(exist(f,'file')==2);
        end

        function exportToExcelFormatWritesEcSheets(testCase)
            % Enzyme-constrained (GECKO) models get two extra export-only
            % sheets, ENZYMES and ENZRXNS, holding the model.ec contents.
            model = testCase.model;
            model.ec.geckoLight = false;
            model.ec.rxns     = model.rxns(1:2);
            model.ec.kcat     = [13.7; 0];
            model.ec.source   = {'brenda'; ''};
            model.ec.notes    = {'note1'; ''};
            model.ec.eccodes  = {'1.1.1.1'; '2.7.1.1;2.7.1.2'};
            model.ec.genes    = model.genes(1:2);
            model.ec.enzymes  = {'P0A1'; 'P0A2'};
            model.ec.mw       = [51000; NaN];
            model.ec.sequence = {'MABC'; 'MDEF'};
            model.ec.concs    = [NaN; 0.5];
            model.ec.rxnEnzMat = [1 2; 0 1];   % R1: P0A1 x1, P0A2 x2; R2: P0A2 x1

            f = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f));
            evalc('exportToExcelFormat(model, f);');

            sheets = sheetnames(f);
            testCase.verifyTrue(any(strcmp(sheets,'ENZYMES')));
            testCase.verifyTrue(any(strcmp(sheets,'ENZRXNS')));

            % NaN mw/conc are written as blanks; the kcat 0 sentinel is kept.
            enz = readcell(f,'Sheet','ENZYMES');
            testCase.verifyEqual(string(enz(1,2:6)),["ID","GENE","MW","SEQUENCE","CONC"]);
            testCase.verifyEqual(string(enz{2,2}),"P0A1");
            testCase.verifyEqual(enz{2,4},51000,'AbsTol',1e-9);

            ecrxn = readcell(f,'Sheet','ENZRXNS');
            testCase.verifyEqual(string(ecrxn(1,2:7)),["ID","KCAT","SOURCE","NOTE","EC-NUMBER","ENZYMES"]);
            testCase.verifyEqual(ecrxn{2,3},13.7,'AbsTol',1e-9);
            testCase.verifyEqual(string(ecrxn{3,6}),"2.7.1.1;2.7.1.2");
            % ENZYMES column: 'enzyme:count' subunit stoichiometry from rxnEnzMat
            testCase.verifyEqual(string(ecrxn{2,7}),"P0A1:1;P0A2:2");
            testCase.verifyEqual(string(ecrxn{3,7}),"P0A2:1");
        end

        function writeReadYAMLRoundTrip(testCase)
            f = [tempname '.yml'];
            testCase.addTeardown(@() delete(f));
            evalc('writeYAMLmodel(testCase.model, f);');
            evalc('m2 = readYAMLmodel(f);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function writeReadYAMLPreservesCompartmentAnnotations(testCase)
            % Compartment annotations must round-trip, and must not be read
            % back as additional compartments.
            m = testCase.model;
            m.compMiriams = cell(numel(m.comps),1);
            m.compMiriams{1}.name  = {'go'};
            m.compMiriams{1}.value = {'GO:0005737'};
            f = [tempname '.yml'];
            testCase.addTeardown(@() delete(f));
            evalc('writeYAMLmodel(m, f);');
            evalc('m2 = readYAMLmodel(f);');
            testCase.verifyEqual(m2.comps, m.comps);
            testCase.verifyEqual(m2.compNames, m.compNames);
            testCase.verifyEqual(m2.compMiriams{1}.name, {'go'});
            testCase.verifyEqual(m2.compMiriams{1}.value, {'GO:0005737'});
        end

        function exportToTabDelimitedWritesFiles(testCase)
            d = [tempname filesep];
            mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            evalc('exportToTabDelimited(testCase.model, d);');
            testCase.verifyTrue(exist(fullfile(d,'excelRxns.txt'),'file')==2);
        end

        function exportForGitWritesRepo(testCase)
            d = [tempname filesep];
            mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            evalc('exportForGit(testCase.model, ''ec'', d, {''yml''});');
            testCase.verifyNotEmpty(dir(fullfile(d,'**','*.yml')));
        end

        function checkFileExistenceFindsFile(testCase)
            f = fullfile(testCase.ravenRoot,'testing','function_tests', ...
                'test_data','ecoli_textbook.mat');
            out = checkFileExistence(f, 1, false, true);
            testCase.verifyNotEmpty(out);
        end

        function cleanSheetTrimsComments(testCase)
            raw = {'#comment', '#comment'; 'a', 'b'; '', ''};
            out = cleanSheet(raw);
            testCase.verifyClass(out, 'cell');
        end

        function writeExcelRoundTrips(testCase)
            % writeExcel generates the .xlsx (Office Open XML) directly, with
            % no external library; read it back with base MATLAB to confirm
            % the sheets, values and types survive.
            f = [tempname '.xlsx'];
            testCase.addTeardown(@() delete(f));
            s(1).name='RXNS'; s(1).header={'ID','VALUE'};
            s(1).data={'r1',1.5; 'r2',-3};
            s(2).name='METS'; s(2).header={'ID'}; s(2).data={'m1'};
            writeExcel(f, s);
            testCase.verifyEqual(sort(string(sheetnames(f))),["METS";"RXNS"]);
            raw = readcell(f,'Sheet','RXNS');
            testCase.verifyEqual(string(raw(1,:)),["ID","VALUE"]);
            testCase.verifyEqual(string(raw{2,1}),"r1");
            testCase.verifyEqual(raw{2,2},1.5,'AbsTol',1e-9);
        end

    end
end
