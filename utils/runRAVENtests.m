function testResults = runRAVENtests
% runRAVENtests
%  Runs all unit tests found at RAVEN/test/unit_tests, and prints output in
%  command window.

ravenPath = findRAVENroot;
curwd = pwd;
cd(fullfile(ravenPath,'testing','unit_tests'));
testResults = runtests(struct2table(dir('*.m')).name);

cd(curwd);