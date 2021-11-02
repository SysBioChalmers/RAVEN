function [success,blastStructure]=testDiamond(suppressWarnings)
% testDiamond
%   Performs a check for DIAMOND functionality in RAVEN. This function is
%   automatically called by checkInstallation function.
%
%   Input:
%   supressWarnings true if warnings should be supressed (opt, default
%                   true)
%
%   Output: 
%   success         true if the test was successful
%   blastStructure	blastStructure resulting from the DIAMOND check
%
%   NOTE: The purpose of the check is to assess whether the
%   homology search can be successfully performed using existing DIAMOND
%   binary. This testing function is completely standalone, only
%   requiring DIAMOND binary and multi-FASTA file from test_data
%   directory. Since one can build DIAMOND database and run homology search
%   using the same binary file (unlike BLAST+), both these functionalities
%   are tested together and it is not possible to run them separately.
%
%   Usage: [success,blastStructure]=testDiamond(suppressWarnings)

if nargin<1
    suppressWarnings=true;
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
else
    dispEM('Unknown OS, exiting.')
    return
end

%Create an empty blastStructure
blastStructure=[];
 
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
copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.fa'),tmpDIR);

%Construct a DIAMOND database
fprintf(['\tdiamond' binEnd ': makedb...\t\t\t\t\t\t\t\t']);
[res, ~]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" makedb --in "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" --db "' tmpDIR '"']);
if res~=0
    fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
    if ~suppressWarnings
        EM=['DIAMOND makedb did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
end
fprintf('OK\n');

%Run a homology search
fprintf(['\tdiamond' binEnd ': blastp...\t\t\t\t\t\t\t\t']);
[res, ~]=system(['"' fullfile(ravenPath,'software','diamond',['diamond' binEnd]) '" blastp --query "' fullfile(tmpDIR,'yeast_galactosidases.fa') '" --out "' outFile '" --db "' tmpDIR '" --more-sensitive --outfmt 6 qseqid sseqid evalue pident length bitscore ppos --threads ' cores ]);
if res~=0
    fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
    if ~suppressWarnings
        EM=['DIAMOND blastp did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
end
fprintf('OK\n');

%Done with the DIAMOND blastp, do the parsing of the text file
blastStructure.fromId='sce';
blastStructure.toId='sco';
A=readtable(outFile,'Delimiter','\t','Format','%s%s%f%f%f%f%f');
blastStructure.fromGenes=A{:,1};
blastStructure.toGenes=A{:,2};
blastStructure.evalue=table2array(A(:,3));
blastStructure.identity=table2array(A(:,4));
blastStructure.aligLen=table2array(A(:,5));
blastStructure.bitscore=table2array(A(:,6));
blastStructure.ppos=table2array(A(:,7));

%Remove the old tempfiles
delete([outFile '*']);
%Remove temporary folder, since homology search is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%If this line is reached then it is assumed that test was successful
success=1;
end
