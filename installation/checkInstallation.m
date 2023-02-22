function checkInstallation(develMode)
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from MATLAB pathlist
%
% Input: 
%   develMode       logical indicating development mode, which includes
%                   testing of binaries that are required to update KEGG
%                   HMMs (opt, default false)
%
%   Usage: checkInstallation(develMode)

if nargin<1
    develMode=false;
end

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

fprintf('\n*** THE RAVEN TOOLBOX ***\n\n');
%Print the RAVEN version if it is not the development version
fprintf([myStr(' > Installing from location',40) '%f']);
fprintf('%s\n',ravenDir)
fprintf([myStr(' > Checking RAVEN release',40) '%f']);
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    currVer = fgetl(fopen(fullfile(ravenDir,'version.txt')));
    fclose('all');
    fprintf([currVer '\n']);
    try
        newVer=strtrim(webread('https://raw.githubusercontent.com/SysBioChalmers/RAVEN/main/version.txt'));
        newVerNum=str2double(strsplit(newVer,'.'));
        currVer=str2double(strsplit(currVer,'.'));
        for i=1:3
            if currVer(i)<newVerNum(i)
                fprintf([myStr('   > Latest RAVEN release available',40) '%f'])
                fprintf(['[\b' newVer,']\b\n'])
                hasGit=exist(fullfile(ravenDir,'.git'),'file');
                if hasGit==7
                    fprintf('     [\bRun git pull in your favourite git client]\b\n')
                    fprintf('     [\bto get the latest RAVEN release]\b\n');
                else
                    fprintf([myStr('     [\bInstructions on how to upgrade]\b',40) '%f'])
                    fprintf('<a href="https://github.com/SysBioChalmers/RAVEN/wiki/Installation#upgrade-to-latest-raven-release">here</a>\n');
                end
                break
            elseif i==3
                fprintf('   > You are running the latest RAVEN release\n');
            end
        end
    catch
        fprintf([myStr('   > Checking for latest RAVEN release',40) '%f'])
        fprintf('[\bFail]\b\n');
        fprintf('     [\bCannot reach GitHub for release info]\b\n');
    end
else
    fprintf('DEVELOPMENT\n');
end

fprintf([myStr(' > Checking MATLAB release',40) '%f'])
fprintf([version('-release') '\n'])
fprintf([myStr(' > Checking system architecture',40) '%f'])
fprintf([computer('arch'),'\n'])

fprintf([myStr(' > Set RAVEN in MATLAB path',40) '%f'])
subpath=regexp(genpath(ravenDir),pathsep,'split'); %List all subdirectories
pathsToKeep=cellfun(@(x) ~contains(x,'.git'),subpath) & cellfun(@(x) ~contains(x,'doc'),subpath);
try
    addpath(strjoin(subpath(pathsToKeep),pathsep));
    fprintf('Pass\n');
    fprintf([myStr(' > Save MATLAB path',40) '%f'])
    try
        savepath
        fprintf('Pass\n')   
    catch
        fprintf('[\bFail]\b\n')
        fprintf(['   You might have to rerun checkInstallation again\n'...
                 '   next time you start up MATLAB\n'])        
    end
catch
    fprintf('[\bFail]\b\n')
end

if isunix
    fprintf([myStr('   > Make binaries executable',40) '%f'])
    status = makeBinaryExecutable();
    if status == 0
        fprintf('Pass\n')
    else
        fprintf('[\bFail]\b\n')
    end
end

%Check if it is possible to parse an Excel file
fprintf('\n=== Model import and export ===\n');
fprintf([myStr(' > Add Java paths for Excel format',40) '%f'])
try
    %Add the required classes to the static Java path if not already added
    addJavaPaths();
    fprintf('Pass\n')
catch
    fprintf('[\bFail]\b\n')
end
fprintf([myStr(' > Checking libSBML version',40) '%f'])
try
    evalc('importModel(fullfile(ravenDir,''tutorial'',''empty.xml''))');
    try
        libSBMLver=OutputSBML; % Only works in libSBML 5.17.0+
        fprintf([libSBMLver.libSBML_version_string '\n']);
    catch
        fprintf('[\bFail]\b\n')
        fprintf('   An older libSBML version was found, update to version 5.17.0 or higher for a significant improvement of model import\n');
    end
catch
    fprintf('[\bFail]\b\n')
    fprintf('   Download libSBML from http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n');
end
fprintf(' > Checking model import and export\n')
[~,res]=evalc("runtests('importExportTests.m');");

fprintf([myStr('   > Import Excel format',40) '%f'])
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > Export Excel format',40) '%f'])
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > Import SBML format',40) '%f'])
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > Export SBML format',40) '%f'])
if res(4).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
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
%     fprintf('[\bFail]\b\n');
% end

fprintf('\n=== Model solvers ===\n');

%Get current solver. Set it to 'none', if it is not set
fprintf(' > Checking for LP solvers\n')
[~,res]=evalc("runtests('solverTests.m');");

fprintf([myStr('   > glpk',40) '%f'])
if res(1).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > gurobi',40) '%f'])
if res(2).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > soplex',40) '%f'])
if res(3).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end

fprintf([myStr('   > cobra',40) '%f'])
if res(4).Passed == 1
    fprintf('Pass\n')
else
    fprintf('[\bFail]\b\n')
end
fprintf([myStr(' > Set RAVEN solver',40) '%f'])
try
    oldSolver=getpref('RAVEN','solver');
    solverIdx=find(strcmp(oldSolver,{'glpk','gurobi','soplex','cobra'}));
catch
    solverIdx=0;
end
% Do not change old solver if functional
if solverIdx~=0 && res(solverIdx).Passed == 1
    fprintf([oldSolver '\n'])
% Order of preference: gurobi > glpk > soplex > cobra
elseif res(2).Passed == 1
    fprintf('gurobi\n')
    setRavenSolver('gurobi');
elseif res(1).Passed == 1
    fprintf('glpk\n')
    setRavenSolver('glpk');
elseif res(3).Passed == 1
    fprintf('soplex\n')
    setRavenSolver('soplex');    
elseif res(4).Passed == 1
    fprintf('cobra\n')
    setRavenSolver('cobra');
else
    fprintf('None, no functional solvers\n')
    fprintf('    The glpk should always be working, check your RAVEN installation to make sure all files are present\n')
end

fprintf('\n=== Essential binary executables ===\n');
fprintf([myStr(' > Checking BLAST+',40) '%f'])
[~,res]=evalc("runtests('blastPlusTests.m');");
res=interpretResults(res);
if res==false
    fprintf('   This is essential to run getBlast()\n')
end

fprintf([myStr(' > Checking DIAMOND',40) '%f'])
[~,res]=evalc("runtests('diamondTests.m');");
res=interpretResults(res);
if res==false
    fprintf('   This is essential to run the getDiamond()\n')
end

fprintf([myStr(' > Checking HMMER',40) '%f'])
[~,res]=evalc("runtests('hmmerTests.m')");
res=interpretResults(res);
if res==false
    fprintf(['   This is essential to run getKEGGModelFromHomology()\n'...
             '   when using a FASTA file as input\n'])
end

if develMode
    fprintf('\n=== Development binary executables ===\n');
    fprintf('NOTE: These binaries are only required when using KEGG FTP dump files in getKEGGModelForOrganism\n');

    fprintf([myStr(' > Checking CD-HIT',40) '%f'])
    [~,res]=evalc("runtests('cdhitTests.m');");
    interpretResults(res);

    fprintf([myStr(' > Checking MAFFT',40) '%f'])
    [~,res]=evalc("runtests('mafftTests.m');");
    interpretResults(res);
end

fprintf('\n=== Compatibility ===\n');
fprintf([myStr(' > Checking function uniqueness',40) '%f'])
checkFunctionUniqueness();

fprintf('\n*** checkInstallation complete ***\n\n');
end

function res = interpretResults(results)
if results.Failed==0 && results.Incomplete==0
    fprintf('Pass\n');
    res=true;
else
    fprintf('[\bFail]\b\n')
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

function status = makeBinaryExecutable()
if ispc
    status = 0; % No need to run on Windows
    return;
end

binDir = fullfile(findRAVENroot(),'software');

binList = {fullfile(binDir,'blast+','blastp');
           fullfile(binDir,'blast+','makeblastdb');
           fullfile(binDir,'cd-hit','cd-hit');
           fullfile(binDir,'diamond','diamond');
           fullfile(binDir,'hmmer','hmmbuild');
           fullfile(binDir,'hmmer','hmmsearch');
           fullfile(binDir,'mafft','mafft-linux64','mafft.bat');
           fullfile(binDir,'mafft','mafft-mac','mafft.bat');};

for i=1:numel(binList)
    [status,cmdout] = system(['chmod +x "' binList{i} '"*']);
    if status ~= 0
        error('Failed to make %s executable: %s ',binList{i},strip(cmdout))
    end
end

end
