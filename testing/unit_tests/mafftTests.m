%run this test case with the command
%results = runtests('mafftTests.m')
function tests = mafftTests
tests = functiontests(localfunctions);
end

function testMafft(testCase)
%This unit test comprises the functionality test for MAFFT in RAVEN:
% 1. Check for resulting file against the expected one.

%%
%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

%Import structure that contains expected MAFFT results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expCdhitMafftOutput.mat'],'expCdhitMafftOutput');

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;

%Run MAFFT multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);

sourceDir = fileparts(which(mfilename));
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

% Define WSL paths
wslPath.fastaFile=getWSLpath([tmpDIR filesep 'yeast_galactosidases.fa']);
wslPath.outFile=getWSLpath(outFile);
wslPath.mafft=getWSLpath(fullfile(ravenPath,'software','mafft','mafft-linux64','mafft.bat'));

%%
%Run protein multi-sequence alignment with MAFFT
if ismac
    [~, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-mac','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' outFile '"']);
elseif isunix
    [~, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-linux64','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' outFile '"']);
elseif ispc
    [~, ~]=system(['wsl "' wslPath.mafft '" --auto --anysymbol --quiet --thread "' num2str(cores) '" --out "' wslPath.outFile '" "' wslPath.fastaFile '"']);
end

%%
%Open actual MAFFT results file
actMafftOutput=importdata(fullfile(outFile));

%Remove the old tempfiles
delete([outFile '*']);

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%%
%Check 1a: Check if files for MAFFT results are the same
verifyEqual(testCase,actMafftOutput,expCdhitMafftOutput);

%Check 1b: Change actual MAFFT results file and check if test fails
actMafftOutput='abc';
verifyNotEqual(testCase,actMafftOutput,expCdhitMafftOutput);
end
