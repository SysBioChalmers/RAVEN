%run this test case with the command
%results = runtests('tinitTests.m')
function tests = tinitTests
    tests = functiontests(localfunctions);
end

function testparsexpTask1List(testCase)
    sourceDir = fileparts(which(mfilename));
    taskStruct = parseTaskList(strcat(sourceDir, '/test_data/test_tasks.txt'));
    taskStructExcel = parseTaskList(strcat(sourceDir, '/test_data/test_tasks.xls'));
    %check that all fields in the first line are what we expect
    
    expTask1.id='ER';
    expTask1.description='Aerobic rephosphorylation of ATP from glucose';
    expTask1.shouldFail=true;
    expTask1.printFluxes=true;
    expTask1.comments='Messed up reaction';
    expTask1.inputs={'O2[s]';'glucose[s]'};
    expTask1.LBin=[23.6;23.6];
    expTask1.UBin=[23.8;23.8];
    expTask1.outputs={'H2O[s]';'CO2[s]'};
    expTask1.LBout=[26.1;26.1];
    expTask1.UBout=[26.2;26.2];
    expTask1.equations={'ATP[c] + H2O[c] => ADP[c] + Pi[c] + H+[c]'};
    expTask1.LBequ=30.2;
    expTask1.UBequ=30.6;
    expTask1.changed={'ATP[a] + H2O[a] => ADP[a] + Pi[a] + H+[a]'};
    expTask1.LBrxn=56.2;
    %check that the check works as expected
    verifyNotEqual(testCase,taskStruct(1),expTask1)
    expTask1.UBrxn=60;%now add the last

    
    verifyEqual(testCase,taskStruct(1),expTask1)
    verifyEqual(testCase,taskStructExcel(1),expTask1)
    
    %check that we have 2 tasks in total
    verifyEqual(testCase,length(taskStruct),2)
    verifyEqual(testCase,length(taskStructExcel),2)

    expTask2.id='BS';
    expTask2.description='ATP de novo synthesis';
    expTask2.shouldFail=false;
    expTask2.printFluxes=false;
    expTask2.comments='';
    expTask2.inputs={'O2[s]';'glucose[s]';'NH3[s]';'Pi[s]'};
    expTask2.LBin=[0;0;0;0];
    expTask2.UBin=[1000;1000;1000;1000];
    expTask2.outputs={'H2O[s]';'CO2[s]';'ATP[c]'};
    expTask2.LBout=[0;0;1];
    expTask2.UBout=[1.1;1.1;1.3];
    expTask2.equations={'ATP[c] + H2O[c] <=> ADP[c] + Pi[c] + H+[c]'};
    expTask2.LBequ=-1000;
    expTask2.UBequ=1000;
    expTask2.changed={};
    expTask2.LBrxn=[];
    expTask2.UBrxn=[];

    verifyEqual(testCase,taskStruct(2),expTask2)
    verifyEqual(testCase,taskStructExcel(2),expTask2)

    %and, check that the two formats produce the same
    verifyEqual(testCase,taskStruct,taskStructExcel)
end