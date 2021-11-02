function [success,blastStructure]=testBlast(testMethod,suppressWarnings)
% testBlast
%   Performs a check for BLAST+ functionality in RAVEN. This function is
%   automatically called by checkInstallation function.
%
%   Input:
%   testMethod      method which should be tested. Available options are
%                   'blastp', 'makeblastdb' and 'both' (opt, default
%                   'both')
%   supressWarnings true if warnings should be supressed (opt, default
%                   true)
%
%   Output:
%   success         true if the test was successful
%   blastStructure	blastStructure resulting from the BLAST+ check. Will be
%                   empty if test was done only for 'makeblastdb'
%
%   NOTE: The purpose of the check is to assess whether the
%   homology search can be successfully performed using existing BLAST+
%   binaries. This testing function is completely standalone, only
%   requiring BLAST+ binaries; multi-FASTA and BLAST database files
%   from test_data directory
%
%   Usage: [success,blastStructure]=testBlast(testMethod,suppressWarnings)

if nargin<1
    testMethod='both';
end
if nargin<2
    suppressWarnings=true;
end

if ~strcmp(testMethod,'blastp') && ~strcmp(testMethod,'makeblastdb') && ~strcmp(testMethod,'both')
    dispEM('testMethod not recognized, exiting.')
    return
end

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

%Create an empty blastStructure
blastStructure=[];

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
copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.fa'),tmpDIR);

if (strcmp(testMethod,'makeblastdb') || strcmp(testMethod,'both'))
    %Construct a BLAST database
    fprintf(['\tmakeblastdb' binEnd '...\t\t\t\t\t\t\t']);
    [res, ~]=system(['"' fullfile(ravenPath,'software','blast+',['makeblastdb' binEnd]) '" -in "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' fullfile(tmpDIR) '" -dbtype prot']);
    if res~=0
        fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
        if ~suppressWarnings
            EM=['makeblastdb did not run successfully, error: ', num2str(res)];
            dispEM(EM,true);
        end
    end
    fprintf('OK\n');
end

if (strcmp(testMethod,'blastp') || strcmp(testMethod,'both'))
    if (strcmp(testMethod,'blastp') && ~strcmp(testMethod,'both'))
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.pdb'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.phr'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.pin'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.pot'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.psq'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.ptf'),tmpDIR);
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.pto'),tmpDIR);
    end
    %Run a homology search
    fprintf(['\tblastp' binEnd '...\t\t\t\t\t\t\t\t']);
    if (strcmp(testMethod,'both'))
        [res, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' outFile '" -db "' fullfile(tmpDIR) '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    else
        [res, ~]=system(['"' fullfile(ravenPath,'software','blast+',['blastp' binEnd]) '" -query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" -out "' outFile '" -db "' fullfile(tmpDIR,'yeast_galactosidases') '" -evalue 10e-5 -outfmt "10 qseqid sseqid evalue pident length bitscore ppos" -num_threads "' cores '"']);
    end
    if res~=0
        fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
        if ~suppressWarnings
            EM=['blastp did not run successfully, error: ', num2str(res)];
            dispEM(EM,true);
        end
    end
    fprintf('OK\n');
    
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

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%If this line is reached then it is assumed that test was successful
success=1;
end
