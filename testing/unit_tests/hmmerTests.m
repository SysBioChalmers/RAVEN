%run this test case with the command
%results = runtests('hmmerTests.m')
function tests = hmmerTests
tests = functiontests(localfunctions);
end

function testHmmer(testCase)
%Smoke test that the bundled hmmsearch binary runs. hmmsearch is used by
%getKEGGModelForOrganism to query a proteome against the prebuilt,
%concatenated KO HMM library in a single search.

%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

if ismac
    binEnd='.mac';
else
    binEnd='';
end
if ispc
    cmd=['wsl "' getWSLpath(fullfile(ravenPath,'software','hmmer','hmmsearch')) '" -h'];
else
    cmd=['"' fullfile(ravenPath,'software','hmmer',['hmmsearch' binEnd]) '" -h'];
end
[status, output]=system(cmd);

%hmmsearch -h prints its usage and exits successfully
verifyEqual(testCase,status,0);
verifyTrue(testCase,contains(output,'hmmsearch'));
end
