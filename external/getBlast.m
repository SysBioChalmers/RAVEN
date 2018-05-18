function blastStructure=getBlast(organismID,fastaFile,modelIDs,refFastaFiles)
% getBlast
%   Performs a bidirectional BLASTp between the organism of interest and a
%   set of template organisms.
%
%   organismID      the id of the organism of interest. This should also
%                   match with the id supplied to getModelFromHomology
%   fastaFile       a FASTA file with the protein sequences for the
%                   organism of interest
%   modelIDs        a cell array of model id:s. These must match the
%                   "model.id" fields in the "models" structure if the output
%                   is to be used with getModelFromHomology
%   refFastaFiles   a cell array with the paths to the corresponding FASTA
%                   files
%   
%   blastStructure  structure containing the bidirectional homology
%                   measurements which are used by getModelFromHomology
%
%   NOTE: This function calls BLASTp to perform a bidirectional homology
%   test between the organism of interest and a set of other organisms
%   using standard settings. If you would like to use other homology
%   measurements, please see getBlastFromExcel.
%
%   Usage: blastStructure=getBlast(organismID,fastaFile,modelIDs,...
%           refFastaFiles)
%
%   Eduard Kerkhoven, 2018-05-18
%

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

% Check that the query and reference fasta files are in the current folder
if isrow(refFastaFiles)
    files=horzcat(refFastaFiles,fastaFile);
else
    files=vertcat(refFastaFiles,fastaFile);
end
for i=1:numel(files)
    if ~(exist(fullfile(cd,files{i}),'file')==2)
        error('FASTA file %s cannot be found',string(files{i}));
    elseif any(strfind(strjoin(files,','),' '))
        error('One or more FASTA files have a space in the filename. Remove this before running getBlast');
    end
end

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
else
    dispEM('Unknown OS, exiting.')
    return
end

% Run BLAST multi-threaded to use all logical cores assigned to MATLAB.
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

[status, ~]=system(['"' fullfile(ravenPath,'software','blast-2.6.0+',['makeblastdb' binEnd]) '" -in "' fastaFile{1} '" -out "' tmpDB '" -dbtype prot']);
if status~=0
    EM=['makeblastdb did not run successfully, error: ', num2str(status)];
    dispEM(EM,true);
end

for i=1:numel(refFastaFiles)
    fprintf(['BLASTing "' modelIDs{i} '" against "' organismID{1} '"..\n']);
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast-2.6.0+',['blastp' binEnd]) '" -query "' refFastaFiles{i} '" -out "' outFile '_' num2str(i) '" -db "' tmpDB '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
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
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast-2.6.0+',['makeblastdb' binEnd]) '" -in "' refFastaFiles{i} '" -out "' tmpDB '" -dbtype prot']);
    if status~=0
        EM=['makeblastdb did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','blast-2.6.0+',['blastp' binEnd]) '" -query "' fastaFile{1} '" -out "' outFile '_r' num2str(i) '" -db "' tmpDB '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
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
        A=importdata([outFile '_' num2str(i)]);
    else
        tempStruct.fromId=organismID{1};
        tempStruct.toId=modelIDs{i-numel(refFastaFiles)};
        A=importdata([outFile '_r' num2str(i-numel(refFastaFiles))]);
    end
    tempStruct.fromGenes=A.textdata(:,1);
    tempStruct.toGenes=A.textdata(:,2);
    tempStruct.evalue=A.data(:,1);
    tempStruct.identity=A.data(:,2);
    tempStruct.aligLen=A.data(:,3);
    tempStruct.bitscore=A.data(:,4);
    tempStruct.ppos=A.data(:,5);
    blastStructure=[blastStructure tempStruct];
end

%Remove the old tempfiles
delete([outFile '*']);
end
