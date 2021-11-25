%run this test case with the command
%results = runtests('mafftTests.m')
function tests = mafftTests
tests = functiontests(localfunctions);
end

function testMafft(testCase)
%This unit test comprises the functionality test for MAFFT in RAVEN:
% 1. MD5 checksum check for MAFFT results file against the expected
%    one.

%%
%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

%Identify the operating system
if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='.exe';
else
    dispEM('Unknown OS, exiting.')
    return
end

%Create empty structures needed for actual MD5 hashes for MAFFT results
actMafftOutputHash={};

expMafftOutputHash='9682b50582529cabd7455c5da943524d';

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

%%
%Run protein multi-sequence alignment with MAFFT
if ismac
    [~, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-mac','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' outFile '"']);
elseif isunix
    [~, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-linux64','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' outFile '"']);
elseif ispc
    [~, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-win','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' outFile '"']);
end

%%
%Calculate MD5 checksum for MAFFT results file
actMafftOutputHash=getMD5Hash(outFile,binEnd);

%Remove the old tempfiles
delete([outFile '*']);

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%%
%Check 1a: Check if MD5 checksums for MAFFT results are the same
verifyEqual(testCase,actMafftOutputHash,expMafftOutputHash);

%Check 1b: Change MD5 checksum and check if test fails
actMafftOutputHash='abc';
verifyNotEqual(testCase,actMafftOutputHash,expMafftOutputHash);
end
