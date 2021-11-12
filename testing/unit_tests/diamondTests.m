%run this test case with the command
%results = runtests('diamondTests.m')
function tests = diamondTests
tests = functiontests(localfunctions);
end

function testDiamond(testCase)
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
else
    dispEM('Unknown OS, exiting.')
    return
end

%Create empty structures needed for actual MD5 hashes and DIAMOND BLAST
%results
actDiamondBlastDbHash={};
actDiamondBlastOutputHash={};
actDiamondBlastStructure=[];

%Create structures that contain expected MD5 hashes and DIAMOND BLAST
%results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expDiamondResults.mat'],'expDiamondBlastOutputHash','expDiamondBlastStructure','expDiamondBlastDbHash');

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;

%Run DIAMOND multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

%Construct a DIAMOND BLAST database
[~, actDiamondBlastDbMsg]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" makedb --in "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" --db "' tmpDIR '"']);

%Run a homology search
[~, ~]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" blastp --query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" --out "' outFile '" --db "' tmpDIR '" --more-sensitive --outfmt 6 qseqid sseqid evalue pident length bitscore ppos --threads ' cores ]);

%Generate actual hashing messages for DIAMOND BLAST database files
switch binEnd
    case '.mac'
        [~, actOutFileHashingMsg]=system(['md5 "' outFile '"']);
    case ''
        [~, actOutFileHashingMsg]=system(['md5sum "' outFile '"']);
    case '.exe'
        [~, actOutFileHashingMsg]=system(['certutil -hashfile "' outFile '" MD5"']);
end

actDiamondBlastDbHash = char(regexp(actDiamondBlastDbMsg,'[a-f0-9]{32}','match'));
actDiamondBlastOutputHash = char(regexp(actOutFileHashingMsg,'[a-f0-9]{32}','match'));

%Do the parsing of the text file
actDiamondBlastStructure.fromId='sce';
actDiamondBlastStructure.toId='sco';
A=readtable(outFile,'Delimiter','\t','Format','%s%s%f%f%f%f%f');
actDiamondBlastStructure.fromGenes=A{:,1};
actDiamondBlastStructure.toGenes=A{:,2};
actDiamondBlastStructure.evalue=table2array(A(:,3));
actDiamondBlastStructure.identity=table2array(A(:,4));
actDiamondBlastStructure.aligLen=table2array(A(:,5));
actDiamondBlastStructure.bitscore=table2array(A(:,6));
actDiamondBlastStructure.ppos=table2array(A(:,7));

%Remove the old tempfiles
delete([outFile '*']);

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);


%Check if MD5 checksums for DIAMOND BLAST database are the same
verifyEqual(testCase,actDiamondBlastDbHash,expDiamondBlastDbHash);

%Change one of the MD5 checksums and check if test fails
actDiamondBlastDbHash=actDiamondBlastOutputHash;
verifyNotEqual(testCase,actDiamondBlastDbHash,expDiamondBlastDbHash);

%Check if MD5 checksums for DIAMOND BLAST results are the same
verifyEqual(testCase,actDiamondBlastOutputHash,expDiamondBlastOutputHash);

%Change MD5 checksum and check if test fails
actDiamondBlastOutputHash=actDiamondBlastStructure;
verifyNotEqual(testCase,actDiamondBlastOutputHash,expDiamondBlastOutputHash);

%Check if DIAMOND BLAST structures are the same
verifyEqual(testCase,actDiamondBlastStructure,expDiamondBlastStructure);

%Modify actual DIAMOND BLAST structure and check if test fails
actDiamondBlastStructure.toId=actDiamondBlastStructure.fromId;
verifyNotEqual(testCase,actDiamondBlastStructure,expDiamondBlastStructure);
end
