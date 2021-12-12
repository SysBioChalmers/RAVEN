%run this test case with the command
%results = runtests('diamondTests.m')
function tests = diamondTests
tests = functiontests(localfunctions);
end

function testDiamond(testCase)
%This unit test comprises several functionality tests for DIAMOND blastp in
%RAVEN:
% 1. MD5 checksum check for DIAMOND database files.
% 2. Non-parsed text check for DIAMOND result files. Although the content
%    of the files is exactly the same, their MD5 hashes are somehow
%    different between the operating systems.
% 3. Check of resulting blastStructure against the expected one. This is
%    done to test DIAMOND blastp results parsing in RAVEN.

%%
%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

%Import structures that contain expected MD5 hashes and DIAMOND blastp
%results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expDiamondResults.mat'],'expBlastStructure','expDiamondReport');

organismID='sce';
fastaFile=fullfile(ravenPath,'test','unit_tests','test_data','yeast_galactosidases.fa');
modelIDs={'hsa' 'afv'};
refFastaFiles={fullfile(ravenPath,'test','unit_tests','test_data','human_galactosidases.fa') fullfile(ravenPath,'test','unit_tests','test_data','aflavus_galactosidases.fa')};

%%
%Run DIAMOND blastp
[actBlastStructure,actDiamondReport]=getDiamond(organismID,fastaFile,modelIDs,refFastaFiles,true,true);

%%
%Test 1a: Check if MD5 checksums for DIAMOND blastp database files are the same
verifyEqual(testCase,actDiamondReport.dbHashes,expDiamondReport.dbHashes);

%Test 1b: Change one of the MD5 checksums and check if test fails
actDiamondReport.dbHashes{1,1}=actDiamondReport.dbHashes{1,2};
verifyNotEqual(testCase,actDiamondReport.dbHashes,expDiamondReport.dbHashes);

%Test 2a: Check if DIAMOND blastp result files are the same
verifyEqual(testCase,actDiamondReport.diamondTxtOutput,expDiamondReport.diamondTxtOutput);

%Test 2b: Change actual DIAMOND blastp result file and check if test fails
actDiamondReport.diamondTxtOutput='empty';
verifyNotEqual(testCase,actDiamondReport.diamondTxtOutput,expDiamondReport.diamondTxtOutput);

%Test 3a: Check if DIAMOND blastp structures are the same
verifyEqual(testCase,actBlastStructure,expBlastStructure);

%Test 3b: Modify actual DIAMOND blastp structure and check if test fails
actBlastStructure(1,1).toId=actBlastStructure(1,1).fromId;
verifyNotEqual(testCase,actBlastStructure,expBlastStructure);
end
