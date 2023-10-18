%run this test case with the command
%results = runtests('miriamTests.m')
function tests = miriamTests
tests = functiontests(localfunctions);
end

function editMiriam_and_extractMiriamTest(testCase)
%This function tests the extractMiriam and editMiriam functions

%%
%Load the expected (i.e. sorted) model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

%Try out all combinations, keep metMiriams field, compare with output
modelTest=editMiriam(model,'met',1,'bigg.metabolite','test','add');
testOut{1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met',1,'bigg.metabolite','test','replace');
testOut{end+1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met',1,'bigg.metabolite','test','fill');
testOut{end+1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met','all','bigg.metabolite','test','add');
testOut{end+1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met',true(1,numel(model.mets)),'bigg.metabolite','test','add');
testOut{end+1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met','3pg_c','bigg.metabolite','test','add');
testOut{end+1}=modelTest.metMiriams;
modelTest=editMiriam(model,'met','3pg_c','bigg.metabolite','3pg','add');
testOut{end+1}=modelTest.metMiriams;

%Load the correct results to compare with
load([sourceDir,'/test_data/miriamTestOutput.mat'])

%Check that the actual model is the same as the expected model
verifyEqual(testCase,testOut,testOutOriginal)
end
