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
% 2. Raw text check for BLAST result files. Although the content of the
%    files are exactly the same, their MD5 hashes are somehow different
%    between the operating systems.
% 3. Check of resulting blastStructure against the expected one. This is
%    done to test BLAST results parsing in RAVEN.

%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='.exe';
    setenv('BLASTDB_LMDB_MAP_SIZE','1000000');
else
    dispEM('Unknown OS, exiting.')
    return
end

%Run BLAST multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create empty structures needed for actual MD5 hashes and BLAST results
actDbHashes.phr={};
actDbHashes.pot={};
actDbHashes.psq={};
actDbHashes.pto={};
actBlastOutput={};
actBlastStructure=[];

%Import structures that contain expected MD5 hashes and BLAST results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expBlastResults.mat'],'expBlastStructure','expDbHashes');
expBlastOutput = importdata([sourceDir,'/test_data/expBlastOutput.txt']);

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;%BLAST results file

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

%Construct a BLAST database
[~, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' fullfile(tmpDIR, 'yeast_galactosidases') '" -dbtype prot']);
%Run a homology search
[~, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' outFile '" -db "' fullfile(tmpDIR,'yeast_galactosidases') '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);

%Generate actual hashing messages for BLAST+ database files
switch binEnd
    case '.mac'
        [~, actPhrHashingMsg]=system(['md5 "' fullfile(tmpDIR,'yeast_galactosidases.phr') '"']);
        [~, actPotHashingMsg]=system(['md5 "' fullfile(tmpDIR,'yeast_galactosidases.pot') '"']);
        [~, actPsqHashingMsg]=system(['md5 "' fullfile(tmpDIR,'yeast_galactosidases.psq') '"']);
        [~, actPtoHashingMsg]=system(['md5 "' fullfile(tmpDIR,'yeast_galactosidases.pto') '"']);
    case ''
        [~, actPhrHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.phr') '"']);
        [~, actPotHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.pot') '"']);
        [~, actPsqHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.psq') '"']);
        [~, actPtoHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.pto') '"']);
    case '.exe'
        [~, actPhrHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.phr') '" MD5"']);
        [~, actPotHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.pot') '" MD5"']);
        [~, actPsqHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.psq') '" MD5"']);
        [~, actPtoHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.pto') '" MD5"']);
end

%Extract actual MD5 hashes for BLAST database files
actDbHashes.phr = char(regexp(actPhrHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.pot = char(regexp(actPotHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.psq = char(regexp(actPsqHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.pto = char(regexp(actPtoHashingMsg,'[a-f0-9]{32}','match'));

%Import BLAST results as raw text
actBlastOutput = importdata(outFile);

%Do the parsing of the BLAST results text file
actBlastStructure.fromId='sce';
actBlastStructure.toId='sco';
A=readtable(outFile,'Delimiter',',','Format','%s%s%f%f%f%f%f');
actBlastStructure.fromGenes=A{:,1};
actBlastStructure.toGenes=A{:,2};
actBlastStructure.evalue=table2array(A(:,3));
actBlastStructure.identity=table2array(A(:,4));
actBlastStructure.aligLen=table2array(A(:,5));
actBlastStructure.bitscore=table2array(A(:,6));
actBlastStructure.ppos=table2array(A(:,7));

%Remove the old tempfiles
delete([outFile '*']);
%Remove temporary folder
[~, ~]=system(['rm "' tmpDIR '" -r']);


%Test 1a: Check if MD5 checksums for BLAST database files are the same
verifyEqual(testCase,actDbHashes,expDbHashes);

%Test 1b: Change one of the MD5 checksums and check if test fails
actDbHashes.phr=actDbHashes.pot;
verifyNotEqual(testCase,actDbHashes,expDbHashes);

%Test 2a: Check if BLAST result files are the same
verifyEqual(testCase,actBlastOutput,expBlastOutput);

%Test 2b: Change actual BLAST result file and check if test fails
actBlastOutput='empty';
verifyNotEqual(testCase,actBlastOutput,expBlastOutput);

%Test 3a: Check if BLAST structures are the same
verifyEqual(testCase,actBlastStructure,expBlastStructure);

%Test 3b: Modify actual BLAST structure and check if test fails
actBlastStructure.toId=actBlastStructure.fromId;
verifyNotEqual(testCase,actBlastStructure,expBlastStructure);
end
