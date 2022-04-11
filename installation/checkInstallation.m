function checkInstallation(develMode)
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from MATLAB pathlist
%
%   Input: 
%   develMode       logical indicating development mode, which includes
%                   testing of binaries that are required to update KEGG
%                   HMMs (opt, default false)
%
%   Usage: checkInstallation(develMode)

if nargin<1
    develMode=false;
end

%Check if RAVEN is in the MATLAB path list
paths=textscan(path,'%s','delimiter', pathsep);
paths=paths{1};

%Get the RAVEN path
ravenDir=findRAVENroot();

fprintf('\n*** THE RAVEN TOOLBOX ***\n\n');
%Print the RAVEN version if it is not the development version
fprintf(' > Checking RAVEN release:\t\t\t\t\t\t\t\t');
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    fprintf([fgetl(fopen(fullfile(ravenDir,'version.txt'))) '\n']);
    fclose('all');
else
    fprintf('DEVELOPMENT\n');
end
fprintf([' > Checking MATLAB release:\t\t\t\t\t\t\t\t' version('-release') '\n']);
fprintf(' > Ensure that RAVEN is on the MATLAB path:\t\t\t\t');
subpath=regexp(genpath(ravenDir),pathsep,'split'); %List all subdirectories
pathsToKeep=cellfun(@(x) isempty(strfind(x,'.git')),subpath) & cellfun(@(x) isempty(strfind(x,'doc')),subpath);
addpath(strjoin(subpath(pathsToKeep),pathsep));
savepath
fprintf('Pass\n');

excelFile=fullfile(ravenDir,'tutorial','empty.xlsx');
xmlFile=fullfile(ravenDir,'tutorial','empty.xml');
matFile=fullfile(ravenDir,'tutorial','empty.mat');
ymlFile=fullfile(ravenDir,'tutorial','empty.yml');

%Get the OS specific binary ending (e.g. exe for Windows)
% binEnd = binaryEnding();

%Check if it is possible to parse an Excel file
fprintf('\n=== Model import and export ===\n');
fprintf(' > Add Java paths for Excel format\t\t\t\t\t\t')
try
    %Add the required classes to the static Java path if not already added
    addJavaPaths();
    fprintf('Pass\n')
catch
    fprintf('Fail\n')
end
fprintf(' > Check libSBML version\t\t\t\t\t\t\t\t')
try
    evalc('importModel(fullfile(ravenDir,''tutorial'',''empty.xml''))');
    try
        libSBMLver=OutputSBML; % Only works in libSBML 5.17.0+
        fprintf([libSBMLver.libSBML_version_string '\n']);
    catch
        fprintf('Fail\n')
        fprintf('   An older libSBML version was found, update to version 5.17.0 or higher for a significant improvement of model import\n');
    end
catch
    fprintf('Fail\n')
    fprintf('   Download libSBML from http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n');
end

fprintf(' > Checking model import and export:\n');

res=runtests('importExportTests.m','OutputDetail',0);
fprintf('   > Import Excel format\t\t\t\t\t\t\t\t')
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
fprintf('   > Export Excel format\t\t\t\t\t\t\t\t')
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
fprintf('   > Import SBML format\t\t\t\t\t\t\t\t\t')
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
fprintf('   > Export SBML format\t\t\t\t\t\t\t\t\t')
if res(4).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

%Check if it is possible to import an SBML model using libSBML

%Check if it is possible to import an YAML model
% fprintf(' > Checking import of model in YAML format:\t\t\t');
% try
%     readYaml(ymlFile,true);
%     fprintf('Pass\n');
% catch
%     fprintf('Fail\n');
% end

fprintf('\n=== Model solvers ===\n');
%Define values for keepSolver and workingSolvers, needed for solver
%functionality check
keepSolver=false;
workingSolvers='';
%Get current solver. Set it to 'none', if it is not set
fprintf(' > Checking for functional solvers:\n');
res=runtests('solverTests.m','OutputDetail',0);
fprintf('   > glpk\t\t\t\t\t\t\t\t\t\t\t\t')
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
fprintf('   > gurobi\t\t\t\t\t\t\t\t\t\t\t\t')
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end
fprintf('   > cobra\t\t\t\t\t\t\t\t\t\t\t\t')
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf(' > Current RAVEN solver preference:\t\t\t\t\t\t')
if ~ispref('RAVEN','solver')
    fprintf('None\n')
else
    oldSolver=getpref('RAVEN','solver');
    fprintf([oldSolver,'\n']);
    solverIdx=strcmp(oldSolver,{'glpk','gurobi','cobra'});
end
fprintf(' > New RAVEN solver preference:\t\t\t\t\t\t\t')
% Order of preference: gurobi > glpk > cobra
if exist('oldSolver','var') && res(solverIdx).Passed == 1
    fprintf([oldSolver ' (unchanged)\n'])
    setRavenSolver(oldSolver);
elseif res(2).Passed == 1
    fprintf('gurobi\n')
    setRavenSolver('gurobi');
elseif res(1).Passed == 1
    fprintf('glpk\n')
    setRavenSolver('glpk');
elseif res(3).Passed == 1
    fprintf('cobra\n')
    setRavenSolver('cobra');
else
    fprintf('None, no functional solvers\n')
    fprintf('    The glpk should always be working, check your RAVEN installation to make sure all files are present\n')
end

fprintf('\n=== Essential binary executables ===\n');
fprintf(' > Checking BLAST+:\t\t\t\t\t\t\t\t\t\t');
res=runtests('blastPlusTests.m','OutputDetail',0);
interpretResults(res);
fprintf(' > Checking DIAMOND:\t\t\t\t\t\t\t\t\t');
res=runtests('diamondTests.m','OutputDetail',0);
interpretResults(res);
fprintf(' > Checking HMMER:\t\t\t\t\t\t\t\t\t\t');
res=runtests('hmmerTests.m','OutputDetail',0);
interpretResults(res);

if develMode
    fprintf('\n=== Development binary executables ===\n');
    fprintf('NOTE: These binaries are only required when using KEGG FTP dump files in getKEGGModelForOrganism\n');
    fprintf(' > Checking CD-HIT:\t\t\t\t\t\t\t\t\t\t');
    res=runtests('cdhitTests.m','OutputDetail',0);
    interpretResults(res);
    fprintf(' > Checking MAFFT:\t\t\t\t\t\t\t\t\t\t');
    res=runtests('mafftTests.m','OutputDetail',0);
    interpretResults(res);
end

fprintf('\n=== Compatibility ===\n');
fprintf(' > Checking uniqueness of RAVEN functions:\t\t\t\t');
checkFunctionUniqueness();

fprintf('\n*** checkInstallation complete ***\n\n');
end

function interpretResults(results)
if results.Failed==0 && results.Incomplete==0
    fprintf('Pass\n');
else
    fprintf('Fail\n')
    fprintf('   Download/compile the binary and rerun checkInstallation\n');
end
end
