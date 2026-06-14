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

function testSBMLExportNestedSubSystems(testCase)
%Regression test: a model where reactions have differing numbers of
%subsystems (nested cell-of-cells, e.g. one reaction in two subsystems and
%others in one) must be exportable to SBML. Previously exportModel failed
%with "Dimensions of arrays being concatenated are not consistent" because
%the subSystems flattening could not concatenate entries of unequal length.
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

%Assign nested subSystems: most reactions get a single subsystem, but a few
%get multiple, which is the scenario that triggered the bug.
nRxns=numel(model.rxns);
model.subSystems=repmat({{'Subsystem A'}},nRxns,1);
model.subSystems{1}={'Subsystem A','Subsystem B'};
model.subSystems{2}={'Subsystem B','Subsystem C'};

tmpFile=fullfile(sourceDir,'testing','unit_tests','test_data','_testNested.xml');
%This call previously errored; verify it completes and round-trips.
evalc('exportModel(model,tmpFile)');
evalc('modelImported=importModel(tmpFile)');
delete(tmpFile);

%Subsystems should be preserved (order within a reaction is not guaranteed)
srt=@(c) sort(reshape(c,1,[]));
for i=1:3
    j=find(strcmp(modelImported.rxns,model.rxns{i}));
    verifyEqual(testCase,srt(modelImported.subSystems{j}),srt(model.subSystems{i}));
end
end

function testYAMLexport(testCase)
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'tutorial','empty.mat'), 'emptyModel');
evalc('writeYAMLmodel(emptyModel,fullfile(sourceDir,''testing'',''unit_tests'',''test_data'',''_test.yml''))');
%File will not be exactly equal as it contains the current date and time,
%so md5 or similar would not work. Just check whether file is reasonably
%sized.
s = dir(fullfile(sourceDir,'testing','unit_tests','test_data','_test.yml'));         
filesize = s.bytes;
verifyTrue(testCase,filesize>1290);
delete(fullfile(sourceDir,'testing','unit_tests','test_data','_test.yml'));
end