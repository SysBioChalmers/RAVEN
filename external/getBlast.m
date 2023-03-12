function [blastStructure,blastReport]=getBlast(organismID,fastaFile,...
    modelIDs,refFastaFiles,developMode,hideVerbose)
% getBlast
%   Performs a bidirectional BLAST between the organism of interest and a
%   set of template organisms
%
%   Input:
%   organismID      the id of the organism of interest. This should also
%                   match with the id supplied to getModelFromHomology
%   fastaFile       a FASTA file with the protein sequences for the
%                   organism of interest
%   modelIDs        a cell array of model ids. These must match the
%                   "model.id" fields in the "models" structure if the
%                   output is to be used with getModelFromHomology
%   refFastaFiles   a cell array with the paths to the corresponding FASTA
%                   files
%   developMode     true if blastReport should be generated that is used
%                   in the unit testing function for BLAST+ (opt, default
%                   false)
%   hideVerbose     true if no status messages should be printed (opt,
%                   default false)
%
%   Output:
%   blastStructure  structure containing the bidirectional homology
%                   measurements that can be used by getModelFromHomology
%   blastReport     structure containing MD5 hashes for FASTA database
%                   files and non-parsed BLAST output data. Will be blank
%                   if developMode is false.
%
%   NOTE: This function calls BLAST+ to perform a bidirectional homology
%   test between the organism of interest and a set of other organisms
%   using standard settings. The only filtering this function does is the
%   removal of hits with an E-value higher than 10e-5. The other homology
%   measurements can be implemented using getBlastFromExcel.
%
%   Usage: [blastStructure,blastReport]=getBlast(organismID,fastaFile,...
%    modelIDs,refFastaFiles,developMode,hideVerbose)

if nargin<5
    developMode = false;
end
if nargin<6
    hideVerbose = false;
end

%Everything should be cell arrays
organismID=convertCharArray(organismID);
fastaFile=convertCharArray(fastaFile);
modelIDs=convertCharArray(modelIDs);
refFastaFiles=convertCharArray(refFastaFiles);

%Create blank structures for results
blastStructure=[];
blastReport.dbHashes.phr={};
blastReport.dbHashes.pot={};
blastReport.dbHashes.psq={};
blastReport.dbHashes.pto={};
blastReport.blastTxtOutput={};

%Get the directory for RAVEN Toolbox
ravenPath=findRAVENroot();

%Generate temporary names for BLAST databases and output files
tmpDB=tempname;
outFile=tempname;

%Check for existence of files. If no full path is specified for a file,
%assume that it is in the current folder
if isrow(refFastaFiles)
    files=horzcat(fastaFile,refFastaFiles);
else
    files=vertcat(fastaFile,refFastaFiles);
end

files=checkFileExistence(files,2); %Copy files to temp dir
fastaFile = files(1);
refFastaFiles = files(2:end);

%Identify the operating system
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

%Create a database for the new organism and blast each of the refFastaFiles
%against it
[status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in ' fastaFile{1} ' -out "' fullfile(tmpDB, 'tmpDB') '" -dbtype prot']);
if developMode
    blastReport.dbHashes.phr{numel(blastReport.dbHashes.phr)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.phr'));
    blastReport.dbHashes.pot{numel(blastReport.dbHashes.pot)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.pot'));
    blastReport.dbHashes.psq{numel(blastReport.dbHashes.psq)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.psq'));
    blastReport.dbHashes.pto{numel(blastReport.dbHashes.pto)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.pto'));
end
if status~=0
    EM=['makeblastdb did not run successfully, error: ', num2str(status)];
    dispEM(EM,true);
end

for i=1:numel(refFastaFiles)
    if ~hideVerbose
        fprintf(['BLASTing "' modelIDs{i} '" against "' organismID{1} '"..\n']);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query ' refFastaFiles{i} ' -out "' outFile '_' num2str(i) '" -db "' fullfile(tmpDB, 'tmpDB') '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    if developMode
        blastReport.blastTxtOutput{numel(blastReport.blastTxtOutput)+1}=importdata([outFile '_' num2str(i)]);
    end
    if status~=0
        EM=['blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
end
delete([tmpDB filesep 'tmpDB*']);

%Then create a database for each of the reference organisms and blast the
%new organism against them
for i=1:numel(refFastaFiles)
    if ~hideVerbose
        fprintf(['BLASTing "' organismID{1} '" against "' modelIDs{i} '"..\n']);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in ' refFastaFiles{i} ' -out "' fullfile(tmpDB, 'tmpDB') '" -dbtype prot']);
    if status~=0
        EM=['makeblastdb did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query ' fastaFile{1} ' -out "' outFile '_r' num2str(i) '" -db "' fullfile(tmpDB, 'tmpDB') '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    if developMode
        blastReport.dbHashes.phr{numel(blastReport.dbHashes.phr)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.phr'));
        blastReport.dbHashes.pot{numel(blastReport.dbHashes.pot)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.pot'));
        blastReport.dbHashes.psq{numel(blastReport.dbHashes.psq)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.psq'));
        blastReport.dbHashes.pto{numel(blastReport.dbHashes.pto)+1}=getMD5Hash(fullfile(tmpDB, 'tmpDB.pto'));
        blastReport.blastTxtOutput{numel(blastReport.blastTxtOutput)+1}=importdata([outFile '_r' num2str(i)]);
    end    
    if status~=0
        EM=['blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    delete([tmpDB filesep 'tmpDB*']);
end

%Done with the BLAST, do the parsing of the text files
for i=1:numel(refFastaFiles)*2
    tempStruct=[];
    if i<=numel(refFastaFiles)
        tempStruct.fromId=modelIDs{i};
        tempStruct.toId=organismID{1};
        A=readtable([outFile '_' num2str(i)],'Delimiter',',','Format','%s%s%f%f%f%f%f');
    else
        tempStruct.fromId=organismID{1};
        tempStruct.toId=modelIDs{i-numel(refFastaFiles)};
        A=readtable([outFile '_r' num2str(i-numel(refFastaFiles))],'Delimiter',',','Format','%s%s%f%f%f%f%f');
    end
    tempStruct.fromGenes=A{:,1};
    tempStruct.toGenes=A{:,2};
    tempStruct.evalue=table2array(A(:,3));
    tempStruct.identity=table2array(A(:,4));
    tempStruct.aligLen=table2array(A(:,5));
    tempStruct.bitscore=table2array(A(:,6));
    tempStruct.ppos=table2array(A(:,7));
    blastStructure=[blastStructure tempStruct];
end

%Remove the old tempfiles
delete([outFile '*']);
%Remove the temp fasta files
delete(files{:})
end
