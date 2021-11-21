%run this test case with the command
%results = runtests('blastPlusTests.m')
function tests = blastPlusTests
tests = functiontests(localfunctions);
end

function testBlastPlus(testCase)
%This unit test comprises several functionality tests for BLAST+ in RAVEN:
% 1. MD5 checksum check for BLAST database files. This check is applied for
%    "phr", "pot", "psq" and "pto" files. The remaining files (i.e. "pdb",
%    "pin" and "ptf") are not compared as these seem to be
%    machine-specific.
% 2. Non-parsed text check for BLAST result files. Although the content of
%    the files is exactly the same, their MD5 hashes are somehow different
%    between the operating systems.
% 3. Check of resulting blastStructure against the expected one. This is
%    done to test BLAST results parsing in RAVEN.

%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

%Import structures that contain expected MD5 hashes and BLAST results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expBlastResults.mat'],'expBlastStructure','expBlastReport');

organismID='sce';
fastaFile=fullfile(ravenPath,'testing','unit_tests','test_data','yeast_galactosidases.fa');
modelIDs={'hsa' 'afv'};
refFastaFiles={fullfile(ravenPath,'testing','unit_tests','test_data','human_galactosidases.fa') fullfile(ravenPath,'testing','unit_tests','test_data','aflavus_galactosidases.fa')};

%Run BLAST
[actBlastStructure,actBlastReport]=getBlast(organismID,fastaFile,modelIDs,refFastaFiles,true,true);


%Test 1a: Check if MD5 checksums for BLAST database files are the same
verifyEqual(testCase,actBlastReport.dbHashes,expBlastReport.dbHashes);

%Test 1b: Change one of the MD5 checksums and check if test fails
actBlastReport.dbHashes.phr{1,1}=actBlastReport.dbHashes.phr{1,2};
verifyNotEqual(testCase,actBlastReport.dbHashes,expBlastReport.dbHashes);

%Test 2a: Check if BLAST result files are the same
verifyEqual(testCase,actBlastReport.blastTxtOutput,expBlastReport.blastTxtOutput);

%Test 2b: Change actual BLAST result file and check if test fails
actBlastReport.blastTxtOutput='empty';
verifyNotEqual(testCase,actBlastReport.blastTxtOutput,expBlastReport.blastTxtOutput);

%Test 3a: Check if BLAST structures are the same
verifyEqual(testCase,actBlastStructure,expBlastStructure);

%Test 3b: Modify actual BLAST structure and check if test fails
actBlastStructure(1,1).toId=actBlastStructure(1,1).fromId;
verifyNotEqual(testCase,actBlastStructure,expBlastStructure);
end
