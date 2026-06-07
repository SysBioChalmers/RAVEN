%run this test case with the command
%results = runtests('importExportTests.m')
function tests = importExportTests
tests = functiontests(localfunctions);
end

function testExcelImport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
excelFile=fullfile(sourceDir,'tutorial','empty.xlsx');
%%Prepare test results, uncomment and run when new reference is needed
% modelExcel=importExcelModel(excelFile);
% sbmlFile=fullfile(sourceDir,'tutorial','empty.xml');
% modelSBML=importModel(sbmlFile);
% save(fullfile(sourceDir,'testing','unit_tests','test_data','importExportResults.mat'),'modelExcel','modelSBML');
evalc('model=importExcelModel(excelFile)'); % Repress warnings
load(fullfile(sourceDir,'testing','unit_tests','test_data','importExportResults.mat'), 'modelExcel');
verifyEqual(testCase,model,modelExcel)
end

function testSBMLImport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
sbmlFile=fullfile(sourceDir,'tutorial','empty.xml');
evalc('model=importModel(sbmlFile)'); % Repress warnings
load(fullfile(sourceDir,'testing','unit_tests','test_data','importExportResults.mat'), 'modelSBML');
verifyEqual(testCase,model,modelSBML)
end

function testYAMLimport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
yamlFile=fullfile(sourceDir,'tutorial','empty.yml');
evalc('model=readYAMLmodel(yamlFile)'); % Repress warnings
load(fullfile(sourceDir,'testing','unit_tests','test_data','importExportResults.mat'), 'modelYAML');
verifyEqual(testCase,model,modelYAML)
end

function testExcelExport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');
evalc('exportToExcelFormat(model,fullfile(sourceDir,''testing'',''unit_tests'',''test_data'',''_test.xlsx''))');
%File will not be exactly equal as it contains the current date and time,
%so md5 or similar would not work. Just check whether file is reasonably
%sized.
s = dir(fullfile(sourceDir,'testing','unit_tests','test_data','_test.xlsx'));         
filesize = s.bytes;
verifyTrue(testCase,filesize>17000);
delete(fullfile(sourceDir,'testing','unit_tests','test_data','_test.xlsx'));
end

function testSBMLExport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');
evalc('exportModel(model,fullfile(sourceDir,''testing'',''unit_tests'',''test_data'',''_test.xml''))');
%File will not be exactly equal as it contains the current date and time,
%so md5 or similar would not work. Just check whether file is reasonably
%sized.
s = dir(fullfile(sourceDir,'testing','unit_tests','test_data','_test.xml'));         
filesize = s.bytes;
verifyTrue(testCase,filesize>18500);
delete(fullfile(sourceDir,'testing','unit_tests','test_data','_test.xml'));
end

function testYAMLexport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
testFile=fullfile(sourceDir,'testing','unit_tests','test_data','_test.yml');
load(fullfile(sourceDir,'tutorial','empty.mat'), 'emptyModel');
evalc('writeYAMLmodel(emptyModel,testFile)');
%File will not be exactly equal as it contains the current date and time,
%so md5 or similar would not work. Sanity-check the size and verify the
%output round-trips back through readYAMLmodel — that's a stronger check
%than a byte-count threshold and isn't tied to quoting / float-format
%decisions in the writer.
s = dir(testFile);
verifyTrue(testCase,s.bytes>1000);
evalc('roundTrip=readYAMLmodel(testFile)');
verifyEqual(testCase,numel(roundTrip.mets),numel(emptyModel.mets));
verifyEqual(testCase,numel(roundTrip.rxns),numel(emptyModel.rxns));
verifyEqual(testCase,roundTrip.id,emptyModel.id);
delete(testFile);
end

function testYAMLpreShimRoundTrip(testCase)
%Regression: the pre-shim ("legacy") RAVEN YAML format — geckoLight
%inside metaData, top-level metabolite SMILES, rxnNotes reaction key,
%integer bounds, double-quoted strings, --- document marker — must
%continue to load via readYAMLmodel even after the writer was aligned
%with cobrapy. The current empty.yml fixture is in the legacy format,
%so we test against it directly.
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
yamlFile=fullfile(sourceDir,'tutorial','empty.yml');
evalc('model=readYAMLmodel(yamlFile)');
verifyTrue(testCase,~isempty(model.id));
verifyTrue(testCase,numel(model.mets)>=4);
verifyTrue(testCase,numel(model.rxns)>=1);
%geckoLight inside metaData populates model.ec when present.
if isfield(model,'ec')
    verifyClass(testCase,model.ec.geckoLight,'logical');
end
end