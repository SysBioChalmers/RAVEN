function [success,blastStructure]=testBlast(fullCheck)
% testBlast
%   Performs a check for BLAST+ functionality in RAVEN. Depending on the
%   parameter settings the user can choose between a quick check for
%   binaries or the thorough testing while building BLAST database and
%   running homology search with BLASTP
%
%   Input:
%   fullCheck       true if the thorough BLAST+ testing should be performed
%                   (opt, default true)
%
%   Output: 
%   success         true if the test was successful, otherwise equal to
%                   zero
%   blastStructure	blastStructure resulting from the thorough BLAST+ check
%
%   NOTE: The purpose of the thorough check is to assess whether the
%   homology search can be successfully performed using existing BLAST+
%   binaries. This testing function is completely standalone, only
%   requiring BLAST+ binaries and multi-FASTA file sce.fa from tutorials
%   directory
%
%   Usage: [success,blastStructure]=testBlast(fullCheck)
%
%   Simonas Marcisauskas, 2019-10-04
%

if nargin<1
    fullCheck=true;
end

%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

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

%Create an empty blastStructure. Even if a quick BLAST+ evaluation is
%considered, blastStructure should still be in the output
blastStructure=[];

if ~fullCheck
    fprintf(['Checking blastp' binEnd '... ']);
    [res,~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '"']);
    if res==1
        fprintf('OK\n');
    else
        fprintf('Not OK! The binary must be recompiled from source before running RAVEN\n');
    end
    fprintf(['Checking makeblastdb' binEnd '... ']);
    [res,~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '"']);
    if res==1
        fprintf('OK\n');
    else
        fprintf('Not OK! The binary must be recompiled from source before running RAVEN\n');
    end
else    
    %Generate temporary names for working directory and outFile
    tmpDB=tempname;
    outFile=tempname;
    
    %Run BLAST multi-threaded to use all logical cores assigned to MATLAB
    cores = evalc('feature(''numcores'')');
    cores = strsplit(cores, 'MATLAB was assigned: ');
    cores = regexp(cores{2},'^\d*','match');
    cores = cores{1};
    
    %Create a temporary folder and copy multi-FASTA file there
    [~, ~]=system(['mkdir "' tmpDB '"']);
    copyfile(fullfile(ravenPath,'tutorial','sce.fa'),tmpDB);
    
    %Construct a BLAST database
    fprintf('Testing makeblastdb... ');
    [res, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' fullfile(tmpDB,'sce.fa') '" -out "' tmpDB '" -dbtype prot']);
    if res~=0
        fprintf('Not OK\n');
        EM=['makeblastdb did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
    fprintf('OK\n');
    
    %Run a homology search
    fprintf('Testing blastp... ');
    [res, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fullfile(tmpDB,'sce.fa') '" -out "' outFile '" -db "' tmpDB '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    if res~=0
        fprintf('Not OK\n');
        EM=['blastp did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
    fprintf('OK\n');
    
    %Remove temporary folder, since homology search is finished
    [~, ~]=system(['rm "' tmpDB '" -r']);
    
    %Done with the BLAST, do the parsing of the text file
    blastStructure.fromId='sce';
    blastStructure.toId='sco';
    A=readtable(outFile,'Delimiter',',','Format','%s%s%f%f%f%f%f');
    blastStructure.fromGenes=A{:,1};
    blastStructure.toGenes=A{:,2};
    blastStructure.evalue=table2array(A(:,3));
    blastStructure.identity=table2array(A(:,4));
    blastStructure.aligLen=table2array(A(:,5));
    blastStructure.bitscore=table2array(A(:,6));
    blastStructure.ppos=table2array(A(:,7));
    
    %Remove the old tempfiles
    delete([outFile '*']);
end

success=1;
end
