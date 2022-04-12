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
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

fprintf('\n*** THE RAVEN TOOLBOX ***\n\n');
%Print the RAVEN version if it is not the development version
fprintf([myStr(' > Checking RAVEN release:',40) '%f'])
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    fprintf([fgetl(fopen(fullfile(ravenDir,'version.txt'))) '\n']);
    fclose('all');
else
    fprintf('DEVELOPMENT\n');
end
fprintf([myStr(' > Checking MATLAB release:',40) '%f'])
fprintf([version('-release') '\n'])
fprintf([myStr(' > Set RAVEN in MATLAB path:',40) '%f'])
subpath=regexp(genpath(ravenDir),pathsep,'split'); %List all subdirectories
pathsToKeep=cellfun(@(x) isempty(strfind(x,'.git')),subpath) & cellfun(@(x) isempty(strfind(x,'doc')),subpath);
addpath(strjoin(subpath(pathsToKeep),pathsep));
savepath
fprintf('Pass\n');

%Check if it is possible to parse an Excel file
fprintf('\n=== Model import and export ===\n');
fprintf([myStr(' > Add Java paths for Excel format:',40) '%f'])
try
    %Add the required classes to the static Java path if not already added
    addJavaPaths();
    fprintf('Pass\n')
catch
    fprintf('Fail\n')
end
fprintf([myStr(' > Check libSBML version:',40) '%f'])
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
fprintf(' > Checking model import and export:\n')
res=runtests('importExportTests.m','OutputDetail',0);

fprintf([myStr('   > Import Excel format:',40) '%f'])
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > Export Excel format:',40) '%f'])
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > Import SBML format:',40) '%f'])
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > Export SBML format:',40) '%f'])
if res(4).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

%Check if it is possible to import an YAML model
% fprintf(' > Checking import of model in YAML format:\t\t\t');
% try
%     readYaml(ymlFile,true);
%     fprintf('Pass\n');
% catch
%     fprintf('Fail\n');
% end

fprintf('\n=== Model solvers ===\n');

%Get current solver. Set it to 'none', if it is not set
fprintf(' > Checking for functional LP solvers:\n')
res=runtests('solverTests.m','OutputDetail',0);

fprintf([myStr('   > glpk:',40) '%f'])
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > gurobi:',40) '%f'])
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > cobra:',40) '%f'])
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf(' > Checking for functional MILP solvers:\n')
res=runtests('fillGapsSmallTests.m','OutputDetail',0);

fprintf([myStr('   > glpk:',40) '%f'])
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > gurobi:',40) '%f'])
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr('   > cobra:',40) '%f'])
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('Fail\n')
end

fprintf([myStr(' > Set RAVEN solver:',40) '%f'])
try
    oldSolver=getpref('RAVEN','solver');
    solverIdx=find(strcmp(oldSolver,{'glpk','gurobi','cobra'}));
catch
    solverIdx=0;
end
% Do not change old solver if functional
if solverIdx~=0 && res(solverIdx).Passed == 1
    fprintf([oldSolver '\n'])
% Order of preference: gurobi > glpk > cobra
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
fprintf([myStr(' > Checking BLAST+:',40) '%f'])
res=runtests('blastPlusTests.m','OutputDetail',0);
interpretResults(res);

fprintf([myStr(' > Checking DIAMOND:',40) '%f'])
res=runtests('diamondTests.m','OutputDetail',0);
interpretResults(res);

fprintf([myStr(' > Checking HMMER:',40) '%f'])
res=runtests('hmmerTests.m','OutputDetail',0);
interpretResults(res);

if develMode
    fprintf('\n=== Development binary executables ===\n');
    fprintf('NOTE: These binaries are only required when using KEGG FTP dump files in getKEGGModelForOrganism\n');

    fprintf([myStr(' > Checking CD-HIT:',40) '%f'])
    res=runtests('cdhitTests.m','OutputDetail',0);
    interpretResults(res);

    fprintf([myStr(' > Checking MAFFT:',40) '%f'])
    res=runtests('mafftTests.m','OutputDetail',0);
    interpretResults(res);
end

fprintf('\n=== Compatibility ===\n');
fprintf([myStr(' > Checking function uniqueness:',40) '%f'])
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

function str = myStr(InputStr,len)
str=InputStr;
lenDiff = len - length(str);
if lenDiff < 0
    warning('String too long');
else
    str = [str blanks(lenDiff)];
end
end
