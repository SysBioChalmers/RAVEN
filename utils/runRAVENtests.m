function testResults = runRAVENtests
% runRAVENtests  Run all RAVEN unit tests.
%
% Runs all unit tests found at RAVEN/test/unit_tests, and prints output in
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
cd(fullfile(ravenPath,'testing','unit_tests'));
testResults = runtests(struct2table(dir('*.m')).name);

cd(curwd);