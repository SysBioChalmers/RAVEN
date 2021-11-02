function success=testMafft(suppressWarnings)
% testCdhit
%   Performs a check for MAFFT functionality in RAVEN. This function is
%   automatically called by checkInstallation function.
%
%   Input:
%   supressWarnings true if warnings should be supressed (opt, default
%                   true)
%
%   Output: 
%   success         true if the test was successful
%
%   NOTE: The purpose of the check is to assess whether the protein
%   clustering can be successfully performed using existing MAFFT binary.
%   This testing function is completely standalone, only requiring MAFFT
%   binary and multi-FASTA file from test_data directory.
%
%   Usage: success=testMafft(suppressWarnings)

if nargin<1
    suppressWarnings=true;
end

%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

%Generate temporary name for working directory
tmpDIR=tempname;

%Run MAFFT multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.fa'),tmpDIR);

%Run MAFFT
fprintf('\tmafft.bat...\t\t\t\t\t\t\t\t');
if ismac
    [res, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-mac','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' fullfile(tmpDIR, 'yeast_galactosidases_msa.fa') '"']);
elseif isunix
    [res, ~]=system(['"' fullfile(ravenPath,'software','mafft','mafft-linux64','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' fullfile(tmpDIR, 'yeast_galactosidases_msa.fa') '"']);
elseif ispc
    [~, message]=system(['"' fullfile(ravenPath,'software','mafft','mafft-win','mafft.bat') '" --auto --anysymbol --thread "' num2str(cores) '" "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" > "' fullfile(tmpDIR, 'yeast_galactosidases_msa.fa') '"']);
    if (contains(message, 'error'))
        res = 1;
    else
        res = 0;
    end
end

if res~=0
    fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
    if ~suppressWarnings
        EM=['mafft did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
end
fprintf('OK\n');

%Remove temporary folder, since homology search is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%If this line is reached then it is assumed that test was successful
success=1;
end
