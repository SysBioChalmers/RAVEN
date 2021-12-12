function results = runTestInOctave(testFunctionName)
% runTestInOctave
%   Runs MATLAB unit_tests in Octave. This function makes a temporary Octave-
%   style version of the test_function, runs the tests, and only gives pass/fail
%   without option of further details.
%    
%   Input:
%   testFunctionName    string with a test function name, as available under
%                       RAVEN/test/unit_tests (e.g. 'blastPlusTests').
%
%   Output:
%   results             structure mimicking a MATLAB TestResult object.
%       Failed          1 if the test failed for whatever reason, otherwise 0.
%       Passed          1 if the test passed, otherwise 0.
%       Incomplete      Always 0.
%       Name            Name of test function and its included test.
%       Details         Empty field.
    
%% Copy test function, without MATLAB's prefixed section with functiontests
tic;
pathOri = which(testFunctionName);
[testDir,~] = fileparts(pathOri);
pathTmp = [testDir filesep 'octaveTempTest.m'];
fidOri = fopen(pathOri,'r');
fidTmp = fopen(pathTmp,'w');
funcCntr = 0;
while ~feof(fidOri)
    tline = fgetl(fidOri);
    if funcCntr < 2
        if startsWithOct(tline,'function')
            funcCntr = funcCntr + 1;
            if funcCntr == 2
                testName=regexprep(tline,'function (.*)\(.*','$1');
                fprintf(fidTmp,'%s\n','function octaveTempTest(testCase)');
            end
        end
    else
        fprintf(fidTmp,'%s\n',tline);
    end
end
fclose('all');

%% Run test (now renamed octaveTempTest)
try
    octaveTempTest(''); % If test throws no errors, it succeeds
    results.Failed=0;
    results.Passed=1;
catch
    warning(['Test ' funcTest ' was not successful']);
    results.Failed=1;
    results.Passed=0;
end
delete(pathTmp);
results.Incomplete=0;
results.Name=[testFunctionName '/' testName];
results.Details=[];
results.Duration=toc;
