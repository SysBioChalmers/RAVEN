function success=testCdhit(suppressWarnings)
% testCdhit
%   Performs a check for CD-HIT functionality in RAVEN. This function is
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
%   clustering can be successfully performed using existing CD-HIT binary.
%   This testing function is completely standalone, only requiring CD-HIT
%   binary and multi-FASTA file from test_data directory.
%
%   Usage: success=testCdhit(suppressWarnings)

if nargin<1
    suppressWarnings=true;
end

success=0;

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

%Generate temporary name for working directory
tmpDIR=tempname;

%Run CD-HIT multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.fa'),tmpDIR);

%Run CD-HIT
fprintf(['\tcd-hit' binEnd '... ']);
[res,~]=system(['"' fullfile(ravenPath,'software','cd-hit',['cd-hit' binEnd]) '" -T "' num2str(cores) '" -i "' fullfile(tmpDIR, 'yeast_galactosidases.fa') '" -o "' fullfile(tmpDIR, 'cdhit_output.fa') '" -c 1.0 -n 5 -M 2000']);
if res~=0
    fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
    if ~suppressWarnings
        EM=['cd-hit did not run successfully, error: ', num2str(res)];
        dispEM(EM,true);
    end
else
    fprintf('OK\n');
end

%Remove temporary folder, since homology search is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%If this line is reached then it is assumed that test was successful
success=1;
end
