function [currVer, installType] = checkInstallation(developMode)
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from MATLAB pathlist
%
% Input: 
%   developMode     logical indicating development mode, which includes
%                   testing of binaries that are required to update KEGG
%                   HMMs (optional, default false). If 'versionOnly' is
%                   specified, only the version is reported as currVer, no
%                   further installation or tests are performed.
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

if nargin<1
    developMode=false;
end
if ischar(developMode) && strcmp(developMode,'versionOnly')
    versionOnly = true;
else
    versionOnly = false;
end

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

installType = 2; % If neither git nor add-on, then ZIP was downloaded
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

if isunix
    fprintf(myStr('   > Make binaries executable',40))
    status = makeBinaryExecutable(ravenDir);
    if status == 0
        fprintf('Pass\n')
    else
        printOrange('Fail\n')
    end
end

%Check if it is possible to parse an Excel file
fprintf('\n=== Model import and export ===\n');
fprintf(myStr(' > Add Java paths for Excel format',40))
try
    %Add the required classes to the static Java path if not already added
    addJavaPaths();
    fprintf('Pass\n')
catch
    printOrange('Fail\n')
end
fprintf(myStr(' > Checking libSBML version',40))
try
    evalc('importModel(fullfile(ravenDir,''tutorial'',''empty.xml''))');
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
fprintf(' > Checking model import and export\n')
[~,res]=evalc("runtests('importExportTests.m');");

fprintf(myStr('   > Import Excel format',40))
if res(1).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
    addList = matlab.addons.installedAddons;
    if any(strcmpi(addList.Name,'Text Analytics Toolbox'))
        fprintf(['   Excel import/export is incompatible with MATLAB Text Analytics Toolbox.\n' ...
                 '   Further instructions => https://github.com/SysBioChalmers/RAVEN/issues/55#issuecomment-1514369299\n'])
    end
end

fprintf(myStr('   > Export Excel format',40))
if res(3).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

fprintf(myStr('   > Import SBML format',40))
if res(2).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

fprintf(myStr('   > Export SBML format',40))
if res(4).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

if res(1).Passed~=1 && res(3).Passed~=1 && exist('vaderSentimentScores.m','file')==2
    fprintf(['   > MATLAB Text Analytics Toolbox found. This should be\n'...
             '     uninstalled if you want to read/write Excel files.\n'...
             '     See RAVEN GitHub Issues page for instructions.\n'])
end

%Check if it is possible to import an YAML model
% fprintf(' > Checking import of model in YAML format:\t\t\t');
% try
%     readYaml(ymlFile,true);
%     fprintf('Pass\n');
% catch
%     printOrange('Fail\n');
% end

fprintf('\n=== Model solvers ===\n');

%Get current solver. Set it to 'none', if it is not set
fprintf(' > Checking for LP solvers\n')
[~,res]=evalc("runtests('solverTests.m');");

fprintf(myStr('   > glpk',40))
if res(1).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

fprintf(myStr('   > gurobi',40))
if res(2).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

fprintf(myStr('   > scip',40))
if res(3).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end

fprintf(myStr('   > cobra',40))
if res(4).Passed == 1
    fprintf('Pass\n')
else
    printOrange('Fail\n')
end
fprintf(myStr(' > Set RAVEN solver',40))
try
    oldSolver=getpref('RAVEN','solver');
    solverIdx=find(strcmp(oldSolver,{'glpk','gurobi','scip','cobra'}));
catch
    solverIdx=0;
end
% Do not change old solver if functional
if solverIdx~=0 && res(solverIdx).Passed == 1
    fprintf([oldSolver '\n'])
% Order of preference: gurobi > glpk > scip > cobra
elseif res(2).Passed == 1
    fprintf('gurobi\n')
    setRavenSolver('gurobi');
elseif res(1).Passed == 1
    fprintf('glpk\n')
    setRavenSolver('glpk');
elseif res(3).Passed == 1
    fprintf('scip\n')
    setRavenSolver('scip');    
elseif res(4).Passed == 1
    fprintf('cobra\n')
    setRavenSolver('cobra');
else
    fprintf('None, no functional solvers\n')
    fprintf('    The glpk should always be working, check your RAVEN installation to make sure all files are present\n')
end

fprintf('\n=== Essential binary executables ===\n');
fprintf(myStr(' > Checking BLAST+',40))
[~,res]=evalc("runtests('blastPlusTests.m');");
res=interpretResults(res);
if res==false
    fprintf('   This is essential to run getBlast()\n')
end

fprintf(myStr(' > Checking DIAMOND',40))
[~,res]=evalc("runtests('diamondTests.m');");
res=interpretResults(res);
if res==false
    fprintf('   This is essential to run the getDiamond()\n')
end

fprintf(myStr(' > Checking HMMER',40))
[~,res]=evalc("runtests('hmmerTests.m')");
res=interpretResults(res);
if res==false
    fprintf(['   This is essential to run getKEGGModelFromHomology()\n'...
             '   when using a FASTA file as input\n'])
end

if developMode
    fprintf('\n=== Development binary executables ===\n');
    fprintf('NOTE: These binaries are only required when using KEGG FTP dump files in getKEGGModelForOrganism\n');

    fprintf(myStr(' > Checking CD-HIT',40))
    [~,res]=evalc("runtests('cdhitTests.m');");
    interpretResults(res);

    fprintf(myStr(' > Checking MAFFT',40))
    [~,res]=evalc("runtests('mafftTests.m');");
    interpretResults(res);
end

fprintf('\n=== Compatibility ===\n');
fprintf(myStr(' > Checking function uniqueness',40))
checkFunctionUniqueness();

fprintf('\n*** checkInstallation complete ***\n\n');
end

function res = interpretResults(results)
if results.Failed==0 && results.Incomplete==0
    fprintf('Pass\n');
    res=true;
else
    printOrange('Fail\n')
    fprintf('   Download/compile the binary and rerun checkInstallation\n');
    res=false;
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

binList = {fullfile(binDir,'blast+','blastp');                      fullfile(binDir,'blast+','blastp.mac');
           fullfile(binDir,'blast+','makeblastdb');                 fullfile(binDir,'blast+','makeblastdb.mac');
           fullfile(binDir,'cd-hit','cd-hit');                      fullfile(binDir,'cd-hit','cd-hit.mac');
           fullfile(binDir,'diamond','diamond');                    fullfile(binDir,'diamond','diamond.mac');
           fullfile(binDir,'hmmer','hmmbuild');                     fullfile(binDir,'hmmer','hmmbuild.mac');
           fullfile(binDir,'hmmer','hmmsearch');                    fullfile(binDir,'hmmer','hmmsearch.mac');
           fullfile(binDir,'GLPKmex','glpkcc.mexa64');              fullfile(binDir,'GLPKmex','glpkcc.mexglx');                 fullfile(binDir,'GLPKmex','glpkcc.mexmaci64');              fullfile(binDir,'GLPKmex','glpkcc.mexmaca64');
           fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexa64'); fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexglx');    fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexmaci64');  fullfile(binDir,'libSBML','TranslateSBML_RAVEN.mexmaca64');
           fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexa64');    fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexglx');       fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexmaci64');     fullfile(binDir,'libSBML','OutputSBML_RAVEN.mexmaca64');
           fullfile(binDir,'mafft','mafft-linux64','mafft.bat');
           fullfile(binDir,'mafft','mafft-mac','mafft.bat');};

for i=1:numel(binList)
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
