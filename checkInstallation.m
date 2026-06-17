function [currVer, installType] = checkInstallation(developMode, checkBinaries)
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from MATLAB pathlist
%
%   NOTE: this function is run before RAVEN has been added to the MATLAB
%   path, so it must not call any other RAVEN functions until it has added
%   RAVEN to the path itself. Its arguments are therefore parsed directly
%   rather than via parseRAVENargs.
%
% Parameters
% ----------
%   developMode     logical indicating development mode, which includes
%                   testing of binaries that are required to update KEGG
%                   HMMs (default false). If 'versionOnly' is
%                   specified, only the version is reported as currVer, no
%                   further installation or tests are performed.
%   checkBinaries   logical whether non-developMode binaries should be
%                   checked for functionality. If false, it overwrites
%                   developMode. Default true.
%
% Output:
%   currVer         current RAVEN version
%   installType     how RAVEN is installed
%                   0:  via git (as .git folder is found)
%                   1:  as MATLAB Add-On
%                   2:  neither of the above, direct download of ZIP file
%                   This matches the installations mentioned in the wiki:
%                   https://github.com/SysBioChalmers/RAVEN/wiki/Installation
%                   0 = advanced / 1 = easy / 2 = medium
%
% Usage: [currVer, installType] = checkInstallation(developMode)

if nargin < 1 || isempty(developMode)
    developMode = false;
end
if nargin < 2 || isempty(checkBinaries)
    checkBinaries = true;
end

if ischar(developMode) && strcmp(developMode,'versionOnly')
    versionOnly = true;
else
    versionOnly = false;
end

%Get the RAVEN path. checkInstallation.m sits in the RAVEN root.
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(ST(I).file);

installType = 2; % If neither git nor add-on, then ZIP was downloaded
addList = matlab.addons.installedAddons;
if isfolder(fullfile(ravenDir,'.git'))
    installType = 0;
elseif any(strcmp(addList.Name,'RAVEN Toolbox'))
    installType = 1;
end

% Do not print first few lines if only version should be reported
if ~versionOnly
    fprintf('\n*** THE RAVEN TOOLBOX ***\n\n');
    %Print the RAVEN version if it is not the development version
    fprintf(myStr(' > Installation type',40));
    switch installType
        case 0
            fprintf('Advanced (via git)\n');
        case 1
            fprintf('Easy (as MATLAB Add-On)\n');
        case 2
            fprintf('Medium (as downloaded ZIP file)\n');
    end
    fprintf(myStr(' > Installing from location',40));
    fprintf('%s\n',ravenDir)
    fprintf(myStr(' > Checking RAVEN release',40));
end
                       
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    currVer = fgetl(fopen(fullfile(ravenDir,'version.txt')));
    fclose('all');
    if ~versionOnly
        fprintf([currVer '\n']);
        try
            newVer=strtrim(webread('https://raw.githubusercontent.com/SysBioChalmers/RAVEN/main/version.txt'));
            newVerNum=str2double(strsplit(newVer,'.'));
            currVerNum=str2double(strsplit(currVer,'.'));
            for i=1:3
                if currVerNum(i)<newVerNum(i)
                    fprintf(myStr('   > Latest RAVEN release available',40))
                    printOrange([newVer,'\n'])
                    switch installType
                        case 0
                            printOrange('     Run git pull in your favourite git client\n')
                            printOrange('     to get the latest RAVEN release\n');
                        case 1
                            printOrange(myStr('     Instructions on how to upgrade',40))
                            fprintf('<a href="https://github.com/SysBioChalmers/RAVEN/wiki/Installation#upgrade-raven-after-easy-installation">here</a>\n');
                        case 2
                            printOrange(myStr('     Instructions on how to upgrade',40))
                            fprintf('<a href="https://github.com/SysBioChalmers/RAVEN/wiki/Installation#upgrade-raven-after-medium-installation">here</a>\n');                            
                    end
                    break
                elseif i==3
                    fprintf('   > You are running the latest RAVEN release\n');
                end
            end
        catch
            fprintf(myStr('   > Checking for latest RAVEN release',40))
            printOrange('Fail\n');
            printOrange('     Cannot reach GitHub for release info\n');
        end
    end
else
    currVer = 'develop';
    if ~versionOnly; fprintf('DEVELOPMENT\n'); end
end
if strcmp(developMode,'versionOnly')
    return;
end

fprintf(myStr(' > Checking MATLAB release',40))
fprintf([version('-release') '\n'])
fprintf(myStr(' > Checking system architecture',40))
fprintf([computer('arch'),'\n'])

fprintf(myStr(' > Set RAVEN in MATLAB path',40))
subpath=regexp(genpath(ravenDir),pathsep,'split'); %List all subdirectories
pathsToKeep=cellfun(@(x) ~contains(x,'.git'),subpath) & cellfun(@(x) ~contains(x,'doc'),subpath);
try
    addpath(strjoin(subpath(pathsToKeep),pathsep));
    fprintf('Pass\n');
    fprintf(myStr(' > Save MATLAB path',40))
    try
        savepath
        fprintf('Pass\n')   
    catch
        printOrange('Fail\n')
        fprintf(['   You might have to rerun checkInstallation again\n'...
                 '   next time you start up MATLAB\n'])        
    end
catch
    printOrange('Fail\n')
end
fprintf(myStr(' > Store RAVEN path as MATLAB pref',40))
try
    setpref('RAVEN','ravenPath',ravenDir);
    fprintf('Pass\n');
catch
    printOrange('Fail\n')
end

if isunix
    fprintf(myStr('   > Make binaries executable',40))
    status = makeBinaryExecutable(ravenDir);
    if status == 0
        fprintf('Pass\n')
    else
        printOrange('Fail\n')
    end
end

%Check the model import and export formats
fprintf('\n=== Model import and export ===\n');
fprintf(myStr(' > Checking libSBML version',40))
model = [];
try
    evalc('model = importModel(fullfile(ravenDir,''tutorial'',''empty.xml''));');
    try
        libSBMLver=OutputSBML_RAVEN; % Only works in libSBML 5.17.0+
        fprintf([libSBMLver.libSBML_version_string '\n']);
    catch
        printOrange('Fail\n')
        fprintf('   An older libSBML version was found, update to version 5.17.0 or higher for a significant improvement of model import\n');
    end
catch
    printOrange('Fail\n')
    fprintf('   Download libSBML from http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n');
end

% Import and export the small "empty" model directly, which is much faster
% than running the full importExportTests test case. A temporary folder
% holds the exported files and is removed afterwards.
fprintf(' > Checking model import and export\n')
tmpDir = tempname; mkdir(tmpDir);

fprintf(myStr('   > Import SBML format',40))
if ~isempty(model)
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

reportCheck('   > Export SBML format', @() exportModel(model, fullfile(tmpDir,'model.xml')));
reportCheck('   > Import YAML format', @() readYAMLmodel(fullfile(ravenDir,'tutorial','empty.yml')));
reportCheck('   > Export YAML format', @() writeYAMLmodel(model, fullfile(tmpDir,'model.yml')));

reportCheck('   > Export Excel format', @() exportToExcelFormat(model, fullfile(tmpDir,'model.xlsx')));

rmdir(tmpDir,'s');

fprintf('\n=== Model solvers ===\n');
fprintf(' > Checking for LP solvers\n')

% Solve a minimal feasible LP with each solver to verify it is functional.
% This is much faster than running the full solverTests test case.
testModel.rxns = {'r1';'r2'};
testModel.mets = {'a'};
testModel.S    = sparse([1 -1]);
testModel.lb   = [0; 0];
testModel.ub   = [1; 1000];
testModel.rev  = [0; 0];
testModel.c    = [0; 1];
testModel.b    = 0;

solvers  = {'glpk','gurobi','scip','cobra'};
solverOK = false(1,numel(solvers));
for i = 1:numel(solvers)
    solverOK(i) = reportCheck(['   > ' solvers{i}], @() testSolver(solvers{i}, testModel));
end

% Keep the current solver if it is still functional, otherwise pick the best
% available one. Order of preference: gurobi > glpk > scip > cobra
fprintf(myStr(' > Set RAVEN solver',40))
oldSolver = '';
try oldSolver = getpref('RAVEN','solver'); catch; end
solverIdx = find(strcmp(oldSolver, solvers));
if ~isempty(solverIdx) && solverOK(solverIdx)
    fprintf([oldSolver '\n'])
    setRavenSolver(oldSolver);
elseif solverOK(2)
    fprintf('gurobi\n')
    setRavenSolver('gurobi');
elseif solverOK(1)
    fprintf('glpk\n')
    setRavenSolver('glpk');
elseif solverOK(3)
    fprintf('scip\n')
    setRavenSolver('scip');
elseif solverOK(4)
    fprintf('cobra\n')
    setRavenSolver('cobra');
else
    fprintf('None, no functional solvers\n')
    fprintf('    The glpk should always be working, check your RAVEN installation to make sure all files are present\n')
end

fprintf('\n=== Essential binary executables ===\n');
if ~checkBinaries
    printOrange('    Skipping check of binary executables\n')
else
    % The BLAST+/DIAMOND/HMMER binaries are downloaded on demand from
    % raven-data (no longer committed to the repository); fetch any that are
    % missing for this platform before checking. This doubles as the up-front
    % "prefetch" for users who want everything ready after installation.
    try
        downloadRavenBinaries({'blast+','diamond','hmmer'});
    catch ME
        printOrange(['    Could not download binaries from raven-data: ' ME.message '\n' ...
            '    (offline? fetch the "*-binaries" RAVEN release for an offline install)\n']);
    end
    % Run each bundled binary with a version/help flag to confirm it
    % executes on this system, rather than running the full test cases.
    if ~reportCheck(' > Checking BLAST+', @() testBinary(ravenDir,'blast+','blastp','-version'))
        fprintf('   This is essential to run getBlast()\n')
    end
    if ~reportCheck(' > Checking DIAMOND', @() testBinary(ravenDir,'diamond','diamond','version'))
        fprintf('   This is essential to run getDiamond()\n')
    end
    if ~reportCheck(' > Checking HMMER', @() testBinary(ravenDir,'hmmer','hmmsearch','-h'))
        fprintf(['   This is essential to run getKEGGModelFromHomology()\n'...
            '   when using a FASTA file as input\n'])
    end
end

fprintf('\n=== Compatibility ===\n');
fprintf(myStr(' > Checking function uniqueness',40))
checkFunctionUniqueness();

fprintf('\n*** checkInstallation complete ***\n');
end

function testSolver(solver, model)
% Set the RAVEN solver and solve a trivial LP; error if it does not return
% an optimal solution.
setRavenSolver(solver);
sol = solveLP(model, 0);
if ~(isfield(sol,'stat') && sol.stat == 1)
    error('Solver %s did not return an optimal solution', solver);
end
end

function testBinary(ravenDir, tool, binName, versionArg)
% Run a bundled binary with versionArg to confirm it executes; error if not.
if ispc
    binEnd = '.exe'; % native Windows builds (raven-data windows-x86_64 ZIPs)
elseif ismac
    binEnd = '.mac';
else
    binEnd = '';
end
cmd = ['"' fullfile(ravenDir,'software',tool,[binName binEnd]) '" ' versionArg];
[status,~] = system(cmd);
if status ~= 0
    error('%s did not execute (exit status %d)', binName, status);
end
end

function ok = reportCheck(label, fcn)
% Print label, run fcn (a function handle) with its output suppressed, and
% report Pass/Fail. Returns true on success.
assert(isa(fcn,'function_handle'));
fprintf(myStr(label,40));
try
    evalc('fcn()');
    ok = true;
    fprintf('Pass\n');
catch
    ok = false;
    printOrange('Fail\n');
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

function status = makeBinaryExecutable(ravenDir)
% This function is required to run when RAVEN is downloaded as MATLAB
% Add-On, in which case the file permissions are not correctly set
if ispc
    status = 0; % No need to run on Windows
    return;
end
binDir = fullfile(ravenDir,'software');

% blast+/diamond/hmmer may be present from an on-demand download
% (downloadRavenBinaries already chmods those) or from the offline bundle;
% chmod whichever unix binaries are present. Missing entries are skipped below.
binList = {fullfile(binDir,'blast+','blastp');                      fullfile(binDir,'blast+','blastp.mac');
           fullfile(binDir,'blast+','makeblastdb');                 fullfile(binDir,'blast+','makeblastdb.mac');
           fullfile(binDir,'diamond','diamond');                    fullfile(binDir,'diamond','diamond.mac');
           fullfile(binDir,'hmmer','hmmsearch');                    fullfile(binDir,'hmmer','hmmsearch.mac');
           fullfile(binDir,'GLPKmex','glpkcc.mexa64');              fullfile(binDir,'GLPKmex','glpkcc.mexglx');                 fullfile(binDir,'GLPKmex','glpkcc.mexmaci64');              fullfile(binDir,'GLPKmex','glpkcc.mexmaca64');
           fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexa64'); fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexglx');    fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexmaci64');  fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexmaca64');
           fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexa64');    fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexglx');       fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexmaci64');     fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexmaca64');};

status = 0;
for i=1:numel(binList)
    if ~isfile(binList{i}); continue; end % skip binaries not present on this platform
    [status,cmdout] = system(['chmod +x "' binList{i} '"']);
    if status ~= 0
        warning('Failed to make %s executable: %s ',binList{i},strip(cmdout))
    end
end
end

function printOrange(stringToPrint)
% printOrange
%   Duplicate of RAVEN/core/printOrange is also kept here, as this function
%   should be able to run before adding RAVEN to the MATLAB path.
try useDesktop = usejava('desktop'); catch, useDesktop = false; end
if useDesktop
    fprintf(['[\b' stringToPrint,']\b'])
else
    fprintf(stringToPrint)
end
end
