function blastStructure=getBlast(organismID,fastaFile,modelIDs,refFastaFiles)
% getBlast
%   Performs a bidirectional BLASTP between the organism of interest and a
%   set of template organisms.
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
%   
%   Output: 
%   blastStructure  structure containing the bidirectional homology
%                   measurements which are used by getModelFromHomology
%
%   NOTE: This function calls BLASTP to perform a bidirectional homology
%   test between the organism of interest and a set of other organisms
%   using standard settings. The only filtering this functions does is the
%   removal of hits with E value higher than 10e-5. If you would like to
%   use other homology measurements, please see getBlastFromExcel.
%
%   Usage: blastStructure=getBlast(organismID,fastaFile,modelIDs,...
%           refFastaFiles)

%Everything should be cell arrays
organismID=cellstr(organismID);
fastaFile=cellstr(fastaFile);
modelIDs=cellstr(modelIDs);
refFastaFiles=cellstr(refFastaFiles);

blastStructure=[];

%Get the directory for RAVEN Toolbox. This may not be the easiest or best
%way to do this
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

%Construct databases and output file
tmpDB=tempname;
outFile=tempname;

%Check for existence of files. If no full path is specified for a file,
%assume that it is in the current folder
if isrow(refFastaFiles)
    files=horzcat(fastaFile,refFastaFiles);
else
    files=vertcat(fastaFile,refFastaFiles);
end

files=checkFileExistence(files,true,false); %No whitespace allowed
fastaFile = files(1);
refFastaFiles = files(2:end);

%Create a database for the new organism and blast each of the refFastaFiles
%against it
if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='';
    setenv('BLASTDB_LMDB_MAP_SIZE','1000000');
else
    dispEM('Unknown OS, exiting.')
    return
end

% Run BLAST multi-threaded to use all logical cores assigned
cores = getNcores();

[status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' fastaFile{1} '" -out "' tmpDB '" -dbtype prot']);
if status~=0
    EM=['makeblastdb did not run successfully, error: ', num2str(status)];
    dispEM(EM,true);
end

for i=1:numel(refFastaFiles)
    fprintf(['BLASTing "' modelIDs{i} '" against "' organismID{1} '"..\n']);
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' refFastaFiles{i} '" -out "' outFile '_' num2str(i) '" -db "' tmpDB '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    if status~=0
        EM=['blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
end
delete([tmpDB '*']);

%Then create a database for each of the reference organisms and blast the
%new organism against them
for i=1:numel(refFastaFiles)
    fprintf(['BLASTing "' organismID{1} '" against "' modelIDs{i} '"..\n']);
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' refFastaFiles{i} '" -out "' tmpDB '" -dbtype prot']);
    if status~=0
        EM=['makeblastdb did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fastaFile{1} '" -out "' outFile '_r' num2str(i) '" -db "' tmpDB '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    delete([tmpDB '*']);
    if status~=0
        EM=['blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
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
end
