%run this test case with the command
%results = runtests('blastPlusTests.m')
function tests = blastPlusTests
tests = functiontests(localfunctions);
end

function testBlastPlus(testCase)
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

%Create empty structures needed for actual MD5 hashes and BLAST results
actDbHashes.pdb={};
actDbHashes.phr={};
actDbHashes.pot={};
actDbHashes.psq={};
actDbHashes.ptf={};
actDbHashes.pto={};
actBlastOutputHash={};
actBlastStructure=[];

%Create structures that contain expected MD5 hashes and BLAST results
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/expBlastResults.mat'],'expBlastOutputHash','expBlastStructure','expDbHashes');

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;

%Run BLAST multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

%Construct a BLAST database
[~, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' fullfile(tmpDIR, 'yeast_galactosidases') '" -dbtype prot']);

%Run a homology search
[~, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' outFile '" -db "' fullfile(tmpDIR,'yeast_galactosidases') '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);

%Generate actual hashing messages for BLAST+ database files

%NOTE: the database file in "pin" format is ignored and not used in
%verifications, because it contains the absolute path to its BLAST+
%database and MD5 checksums are therefore always different between
%actual and expected. It is therefore assumed that MD5 checksum
%comparison for the other 6 files comprising BLAST+ database is enough
%to test makeblastdb functionality.
switch binEnd
    case '.mac'
        [~, actPdbHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.pdb') '"']);
        [~, actPhrHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.phr') '"']);
        [~, actPotHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.pot') '"']);
        [~, actPsqHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.psq') '"']);
        [~, actPtfHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.ptf') '"']);
        [~, actPtoHashingMsg]=system(['md5 -s "' fullfile(tmpDIR,'yeast_galactosidases.pto') '"']);
        [~, actOutFileHashingMsg]=system(['md5 -s "' outFile '"']);
    case ''
        [~, actPdbHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.pdb') '"']);
        [~, actPhrHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.phr') '"']);
        [~, actPotHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.pot') '"']);
        [~, actPsqHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.psq') '"']);
        [~, actPtfHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.ptf') '"']);
        [~, actPtoHashingMsg]=system(['md5sum "' fullfile(tmpDIR,'yeast_galactosidases.pto') '"']);
        [~, actOutFileHashingMsg]=system(['md5sum "' outFile '"']);
    case '.exe'
        [~, actPdbHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.pdb') '" MD5"']);
        [~, actPhrHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.phr') '" MD5"']);
        [~, actPotHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.pot') '" MD5"']);
        [~, actPsqHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.psq') '" MD5"']);
        [~, actPtfHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.ptf') '" MD5"']);
        [~, actPtoHashingMsg]=system(['certutil -hashfile "' fullfile(tmpDIR,'yeast_galactosidases.pto') '" MD5"']);
        [~, actOutFileHashingMsg]=system(['certutil -hashfile "' outFile '" MD5"']);
end

actDbHashes.pdb = char(regexp(actPdbHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.phr = char(regexp(actPhrHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.pot = char(regexp(actPotHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.psq = char(regexp(actPsqHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.ptf = char(regexp(actPtfHashingMsg,'[a-f0-9]{32}','match'));
actDbHashes.pto = char(regexp(actPtoHashingMsg,'[a-f0-9]{32}','match'));
actBlastOutputHash = char(regexp(actOutFileHashingMsg,'[a-f0-9]{32}','match'));

%Do the parsing of the text file
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

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);


%Check if MD5 checksums for BLAST database are the same
verifyEqual(testCase,actDbHashes,expDbHashes);

%Change one of the MD5 checksum and check if test fails
actDbHashes.pdb=actDbHashes.phr;
verifyNotEqual(testCase,actDbHashes,expDbHashes);

%Check if MD5 checksums for BLAST results are the same
verifyEqual(testCase,actBlastOutputHash,expBlastOutputHash);

%Change MD5 checksum and check if test fails
actBlastOutputHash=actDbHashes.phr;
verifyNotEqual(testCase,actBlastOutputHash,expBlastOutputHash);

%Check if BLAST structures are the same
verifyEqual(testCase,actBlastStructure,expBlastStructure);

%Modify actual BLAST structure and check if test fails
actBlastStructure.toId=actBlastStructure.fromId;
verifyNotEqual(testCase,actBlastStructure,expBlastStructure);
end
