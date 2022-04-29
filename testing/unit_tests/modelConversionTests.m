%run this test case with the command
%results = runtests('modelConversionTests.m')
function tests = modelConversionTests
tests = functiontests(localfunctions);
end

function ravenCobraWrapperTest(testCase)
%This function tests the ravenCobraWrapper function by going through a full
%cycle RAVEN -> COBRA -> RAVEN, and the model should be (near) identical
%(only annotation field should be lost).

%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

%Prevent writing of messages to command window
evalc('modelCobra=ravenCobraWrapper(model);')
evalc('modelRaven=ravenCobraWrapper(modelCobra);')

%We know that the annotation field is lost
modelNoAnnot=rmfield(model,'annotation');

%Check that the actual model is the same as the expected model
verifyEqual(testCase,modelRaven,modelNoAnnot)
end
