function [blastStructure,diamondReport]=getDiamond(organismID,fastaFile,...
    modelIDs,refFastaFiles,develMode,hideVerbose)
% getDiamond
%   Uses DIAMOND to perform a bidirectional BLAST between the organism
%   of interest and a set of template organisms
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
%   develMode       true if blastReport should be generated that is used
%                   in the unit testing function for DIAMOND (opt, default
%                   false)
%   hideVerbose     true if no status messages should be printed (opt,
%                   default false)
%
%   Output:
%   blastStructure  structure containing the bidirectional homology
%                   measurements which are used by getModelFromHomology
%   diamondReport   structure containing MD5 hashes for FASTA database
%                   files and non-parsed BLAST output data. Will be blank
%                   if develMode is false.
%
%   NOTE: This function calls DIAMOND to perform a bidirectional homology
%   search between the organism of interest and a set of other organisms
%   using the '--more-sensitive' setting from DIAMOND. For the most
%   sensitive results, the use of getBlast() is adviced, however,
%   getDiamond() is a fast alternative (>15x faster). The blastStructure
%   generated is in the same format as those obtained from getBlast().
%
%   Usage: [blastStructure,diamondReport]=getDiamond(organismID,fastaFile,...
%    modelIDs,refFastaFiles,develMode,hideVerbose)

if nargin<5
    develMode = false;
end
if nargin<6
    hideVerbose = false;
end

%Everything should be cell arrays
organismID=cellstr(organismID);
fastaFile=cellstr(fastaFile);
modelIDs=cellstr(modelIDs);
refFastaFiles=cellstr(refFastaFiles);

%Create blank structures for results
blastStructure=[];
diamondReport.dbHashes={};
diamondReport.diamondTxtOutput={};

%Get the directory for RAVEN Toolbox. This may not be the easiest or best
%way to do this
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

%Generate temporary names for DIAMOND databases and output files
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
else
    dispEM('Unknown OS, exiting.')
    return
end

cores = getNcores();
% Run DIAMOND multi-threaded to use all logical cores assigned

%Create a database for the new organism and blast each of the refFastaFiles
%against it
[status, message]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" makedb --in "' fastaFile{1} '" --db "' fullfile(tmpDB) '"']);
if develMode
    diamondReport.dbHashes{numel(diamondReport.dbHashes)+1} = char(regexp(message,'[a-f0-9]{32}','match'));
end
if status~=0
    EM=['DIAMOND makedb did not run successfully, error: ', num2str(status)];
    dispEM(EM,true);
end

for i=1:numel(refFastaFiles)
    if ~hideVerbose
        fprintf(['Running DIAMOND blastp with "' modelIDs{i} '" against "' organismID{1} '"..\n']);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" blastp --query "' refFastaFiles{i} '" --out "' outFile '_' num2str(i) '" --db "' fullfile(tmpDB) '" --more-sensitive --outfmt 6 qseqid sseqid evalue pident length bitscore ppos --threads ' cores ]);
    if develMode
        diamondReport.diamondTxtOutput{numel(diamondReport.diamondTxtOutput)+1}=importdata([outFile '_' num2str(i)]);
    end
    if status~=0
        EM=['DIAMOND blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
end
delete([tmpDB filesep 'tmpDB*']);

%Then create a database for each of the reference organisms and blast the
%new organism against them
for i=1:numel(refFastaFiles)
    if ~hideVerbose
        fprintf(['Running DIAMOND blastp with "' organismID{1} '" against "' modelIDs{i} '"..\n']);
    end
    [status, message]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" makedb --in "' refFastaFiles{i} '" --db "' fullfile(tmpDB) '"']);
    if status~=0
        EM=['DIAMOND makedb did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    [status, ~]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" blastp --query "' fastaFile{1} '" --out "' outFile '_r' num2str(i) '" --db "' fullfile(tmpDB) '" --more-sensitive --outfmt 6 qseqid sseqid evalue pident length bitscore ppos --threads ' cores]);
    if develMode
        diamondReport.dbHashes{numel(diamondReport.dbHashes)+1} = char(regexp(message,'[a-f0-9]{32}','match'));
        diamondReport.diamondTxtOutput{numel(diamondReport.diamondTxtOutput)+1}=importdata([outFile '_r' num2str(i)]);
    end
    if status~=0
        EM=['DIAMOND blastp did not run successfully, error: ', num2str(status)];
        dispEM(EM,true);
    end
    delete([tmpDB filesep 'tmpDB*']);
end

%Done with the DIAMOND blastp, do the parsing of the text files
for i=1:numel(refFastaFiles)*2
    tempStruct=[];
    if i<=numel(refFastaFiles)
        tempStruct.fromId=modelIDs{i};
        tempStruct.toId=organismID{1};
        A=readtable([outFile '_' num2str(i)],'Delimiter','\t','Format','%s%s%f%f%f%f%f');
    else
        tempStruct.fromId=organismID{1};
        tempStruct.toId=modelIDs{i-numel(refFastaFiles)};
        A=readtable([outFile '_r' num2str(i-numel(refFastaFiles))],'Delimiter','\t','Format','%s%s%f%f%f%f%f');
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
