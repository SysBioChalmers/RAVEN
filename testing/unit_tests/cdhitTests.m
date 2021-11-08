%run this test case with the command
%results = runtests('cdhitTests.m')
function tests = cdhitTests
tests = functiontests(localfunctions);
end

function testCdhit(testCase)
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

%Create empty structures needed for actual MD5 hashes for CD-HIT results
actCdhitOutputHash={};

expCdhitOutputHash='9682b50582529cabd7455c5da943524d';

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;

%Run CD-HIT multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);

sourceDir = fileparts(which(mfilename));
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

%Run protein clustering with CD-HIT
[~, ~]=system(['"' fullfile(ravenPath,'software','cd-hit',['cd-hit' binEnd]) '" -T "' num2str(cores) '" -i "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" -o "' outFile '" -c 1.0 -n 5 -M 2000']);

%Generate actual hashing messages for CD-HIT results
switch binEnd
    case '.mac'
        [~, actOutFileHashingMsg]=system(['md5 -s "' outFile '"']);
    case ''
        [~, actOutFileHashingMsg]=system(['md5sum "' outFile '"']);
    case '.exe'
        [~, actOutFileHashingMsg]=system(['certutil -hashfile "' outFile '" MD5"']);
end

actCdhitOutputHash = char(regexp(actOutFileHashingMsg,'[a-f0-9]{32}','match'));

%Remove the old tempfiles
delete([outFile '*']);

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);


%Check if MD5 checksums for CD-HIT results are the same
verifyEqual(testCase,actCdhitOutputHash,expCdhitOutputHash);

%Change MD5 checksum and check if test fails
actCdhitOutputHash='abc';
verifyNotEqual(testCase,actCdhitOutputHash,expCdhitOutputHash);
end
