%run this test case with the command
%results = runtests('modelAbilitiesTests.m')
function tests = modelAbilitiesTests
tests = functiontests(localfunctions);
end

function canConsumeTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

testOut=canConsume(model,{'r5p_c';'ru5p_DASH_D_c';'s7p_c';'succ_c';'succ_e';'succoa_c';'xu5p_DASH_D_c'});

testCheck=[true;true;true;true;true;false;true];

verifyEqual(testCase,testOut,testCheck)
end

function canProduceTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

testOut=canProduce(model,{'r5p_c';'ru5p_DASH_D_c';'s7p_c';'succ_c';'succ_e';'succoa_c';'xu5p_DASH_D_c'});

testCheck=[true;true;true;true;true;false;true];

verifyEqual(testCase,testOut,testCheck)
end

function checkProductionTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

[testOut.np, testOut.npn, testOut.nfpm,testOut.mtc,~]=checkProduction(model,true,model.comps,false);

testCheck.np=[10;13;16;17;21;30;32;37;49;50;51;52;53;64;65;71];
testCheck.npn={'Acetyl-CoA[c]';'ADP[c]';'AMP[c]';'ATP[c]';'Coenzyme A[c]';'D-Fructose[e]';'Fumarate[e]';'L-Glutamine[e]';'L-Malate[e]';'NAD[c]';'NADH[c]';'NADP[c]';'NADPH[c]';'Ubiquinone-8[c]';'Ubiquinol-8[c]';'Succinyl-CoA[c]'};
testCheck.nfpm=[true,false,false,false,true,false,false,false,false,false,false,false,false,false,false,true;false,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false;false,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false;false,true,true,true,false,false,false,false,false,false,false,false,false,false,false,false;true,false,false,false,true,false,false,false,false,false,false,false,false,false,false,true;false,false,false,false,false,true,false,false,false,false,false,false,false,false,false,false;false,false,false,false,false,false,true,false,false,false,false,false,false,false,false,false;false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false;false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false;false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false;false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false;false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false;false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false;false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false;false,false,false,false,false,false,false,false,false,false,false,false,false,true,true,false;true,false,false,false,true,false,false,false,false,false,false,false,false,false,false,true];
testCheck.mtc={'Acetyl-CoA[c] (connects 3 metabolites)';'ADP[c] (connects 3 metabolites)';'NAD[c] (connects 2 metabolites)';'NADP[c] (connects 2 metabolites)';'Ubiquinone-8[c] (connects 2 metabolites)';'D-Fructose[e] (connects 1 metabolites)';'Fumarate[e] (connects 1 metabolites)';'L-Glutamine[e] (connects 1 metabolites)';'L-Malate[e] (connects 1 metabolites)'};

verifyEqual(testCase,testOut,testCheck)
end

function consumeSomethingTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

[testOut.s, testOut.m]=consumeSomething(model,{'r5p_c';'ru5p_DASH_D_c';'s7p_c';'succ_c';'succ_e';'succoa_c';'xu5p_DASH_D_c'});
% Simplify output for later check, no need for high precision
% testOut.s(abs(testOut.s)<1e-10)=0;
% testOut.s=round(testOut.s,4)
testCheck.s=[0;0;1.4591;0;0;1.4591;0;0;0;0;8.39;4.0126;0;1.4591;0;2.9183;0;1.4591;0;1.4591;0;0;1.4591;0;0;0;0;0;0;0;1.4591;0;0;0;0;-1.4591;1.4591;0;0;0;0;0;0;0;0;0;0;0;-1.4591;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;2.9183;0;0;-1.4591;1.4591;0;0;0;-1.4591;0;1.4591;1.4591;0;0;0;-1.4591;1.4591;0;0;0;0;0;0;0;0;0;0;0;0];
testCheck.m=[33];

verifyEqual(testCase,testOut,testCheck,'AbsTol',1e-4)
end

function makeSomethingTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

[testOut.s, testOut.m]=makeSomething(model,{'r5p_c';'ru5p_DASH_D_c';'s7p_c';'succ_c';'succ_e';'succoa_c';'xu5p_DASH_D_c'});
% Simplify output for later check, no need for high precision
% testOut.s(abs(testOut.s)<1e-10)=0;
% testOut.s=round(testOut.s,4);

testCheck.s=[0;0;1.8644;0;0;0;0;0;0;0;8.39;4.6611;0;0;0;3.7289;0;1.8644;0;0;0;0;0;0;0;0;0;-0.93220;0;0;0;0;0;0;0;-1.8644;0;0;0;0.93220;0;0;0;0;0;0;0;0;-1.8644;0.93220;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;3.7289;0;0;-1.8644;1.8644;0.93220;0;-0.93220;-1.8644;0;1.8644;0;0;0;0;-1.8644;0.93220;0;0;0;0;0;0;0;0;0;0;0;0.93220];
testCheck.m=[6];

verifyEqual(testCase,testOut,testCheck,'AbsTol',1e-4)
end

function getElementalBalanceTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

testOut=getElementalBalance(model,1);

testCheck.balanceStatus=1;
testCheck.elements.abbrevs={'C';'N';'O';'S';'P';'H'};
testCheck.elements.names={'carbon';'nitrogen';'oxygen';'sulfur';'phosphorus';'hydrogen'};
testCheck.leftComp=[44,14,31,1,5,62];
testCheck.rightComp=[44,14,31,1,5,62];

verifyEqual(testCase,testOut,testCheck)
end

function getAllowedBoundsTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');


[testOut.min, testOut.max, testOut.flag]=getAllowedBounds(model,[1:10]);
% Simplify output for later check, no need for high precision
% testOut.min(abs(testOut.min)<1e-10)=0;
% testOut.min=single(testOut.min);
% testOut.max(abs(testOut.max)<1e-10)=0;
% testOut.max=single(testOut.max);

testCheck.min=[-20;0;0;0;0;0;-166.610;0;0;0];
testCheck.max=[0;20;20;20;20;20;0;20;10;20];
testCheck.flag=[1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1];

verifyEqual(testCase,testOut,testCheck,'AbsTol',1e-4)
end

function getEssentialRxnsTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

[testOut.r, testOut.i]=getEssentialRxns(model);

testCheck.r={'EX_glc';'GLCpts'};
testCheck.i=[28;50];
verifyEqual(testCase,testOut,testCheck)
end
