function testResults = runRAVENtests
% runRAVENtests  Run all RAVEN function tests.
%
% Runs all unit tests found at RAVEN/test/function_tests, and prints output in
% the command window.
%
% Returns
% -------
% testResults : matlab.unittest.TestResult
%     array of test results returned by runtests.
%
% Examples
% --------
%     testResults = runRAVENtests;

ravenPath = findRAVENroot;
curwd = pwd;
cd(fullfile(ravenPath,'testing','function_tests'));

% Collect the concrete test classes only. RavenTestCase is an abstract base
% class that the test classes subclass, and runtests cannot build a suite
% from an abstract class, so it must be excluded from the file list.
testFiles = struct2table(dir('*.m')).name;
testFiles(strcmp(testFiles,'RavenTestCase.m')) = [];
testResults = runtests(testFiles);

cd(curwd);