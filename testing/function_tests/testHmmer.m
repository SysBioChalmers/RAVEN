function success=testHmmer(testMethod,suppressWarnings)
% testHmmer
%   Performs a check for HMMER functionality in RAVEN. This function is
%   automatically called by checkInstallation function.
%
%   Input:
%   testMethod      method which should be tested. Available options are
%                   'hmmsearch', 'hmmbuild' and 'both' (opt, default
%                   'both')
%   supressWarnings true if warnings should be supressed (opt, default
%                   true)
%
%   Output:
%   success         true if the test was successful
%
%   NOTE: The purpose of the check is to assess whether the
%   homology search can be successfully performed using existing HMMER
%   binaries. This testing function is completely standalone, only
%   requiring HMMER binaries; multi-FASTA and HMM files from test_data
%   directory
%
%   Usage: success=testHmmer(testMethod,suppressWarnings)

if nargin<1
    testMethod='both';
end
if nargin<2
    suppressWarnings=true;
end

if ~strcmp(testMethod,'hmmsearch') && ~strcmp(testMethod,'hmmbuild') && ~strcmp(testMethod,'both')
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
else
    dispEM('Unknown OS, exiting.')
    return
end

%Generate temporary name for working directory
tmpDIR=tempname;

%Run BLAST multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);
copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.fa'),tmpDIR);

if (strcmp(testMethod,'hmmbuild') || strcmp(testMethod,'both'))
    %Train a hidden Markov model
    fprintf(['\thmmbuild' binEnd '...\t\t\t\t\t\t\t']);
    [res,~]=system(['"' fullfile(ravenPath,'software','hmmer',['hmmbuild' binEnd]) '" --cpu "' num2str(cores) '" "' fullfile(tmpDIR,'yeast_galactosidases.hmm') '" "' fullfile(tmpDIR,'yeast_galactosidases.fa') '"']);
    if res~=0
        fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
        if ~suppressWarnings
            EM=['hmmbuild did not run successfully, error: ', num2str(res)];
            dispEM(EM,true);
        end
    end
    fprintf('OK\n');
end

if (strcmp(testMethod,'hmmsearch') || strcmp(testMethod,'both'))
    if (strcmp(testMethod,'hmmsearch') && ~strcmp(testMethod,'both'))
        copyfile(fullfile(ravenPath,'testing','function_tests','test_data','yeast_galactosidases.hmm'),tmpDIR);
    end
    %Run a homology search
    fprintf(['\thmmsearch' binEnd '...\t\t\t\t\t\t\t\t']);
    [res, ~]=system(['"' fullfile(ravenPath,'software','hmmer',['hmmsearch' binEnd]) '" --cpu "' num2str(cores) '" "' fullfile(tmpDIR,'yeast_galactosidases.hmm') '" "' fullfile(tmpDIR,'yeast_galactosidases.fa') '"']);
    if res~=0
        fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
        if ~suppressWarnings
            EM=['hmmsearch did not run successfully, error: ', num2str(res)];
            dispEM(EM,true);
        end
    end
    fprintf('OK\n');
end

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);

%If this line is reached then it is assumed that test was successful
success=1;
end
