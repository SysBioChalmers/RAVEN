function checkInstallation()
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from MATLAB pathlist
%
%   Usage: checkInstallation()

%Check if RAVEN is in the path list
paths=textscan(path,'%s','delimiter', pathsep);
paths=paths{1};

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(ST(I).file);

currentWd = pwd; % Path not necessarily added yet, navigate to location.
cd([ravenDir filesep 'src' filesep 'utilities' filesep 'octavematlab'])
if isOctave
    mORo = 'Octave';
else
    mORo = 'MATLAB';
end
cd(currentWd);

%Print the RAVEN version if it is not the development version
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    fprintf(['\n*** THE RAVEN TOOLBOX v.' fgetl(fopen(fullfile(ravenDir,'version.txt'))) ' ***\n\n']);
else
    fprintf('\n*** THE RAVEN TOOLBOX - DEVELOPMENT VERSION ***\n\n');
end
fprintf(' > Checking environment...  %s.\n',mORo);
fprintf(' > Checking if RAVEN is on the %s path...  ',mORo);
if ismember(ravenDir,paths)
    fprintf('OK.\n');
else
    subpath=regexp(genpath(ravenDir),pathsep,'split'); %List all subdirectories
    pathsToKeep=cellfun(@(x) isempty(strfind(x,'.git')),subpath) & cellfun(@(x) isempty(strfind(x,'doc')),subpath);
    addpath(strjoin(subpath(pathsToKeep),pathsep));
    savepath
    fprintf('OK, just added.\n');
end

excelFile=fullfile(ravenDir,'tutorial','empty.xlsx');
xmlFile=fullfile(ravenDir,'tutorial','empty.xml');
matFile=fullfile(ravenDir,'tutorial','empty.mat');
ymlFile=fullfile(ravenDir,'tutorial','empty.yml');

%Check if it is possible to parse an Excel file
fprintf('\n=== Model import and export ===\n');
fprintf(' > Checking export of model in Excel format...  ');
if isOctave
    fprintf('Skipped, MATLAB only.\n')
else
    try
        %Add the required classes to the static Java path if not already added
        %addJavaPaths();
        importExcelModel(excelFile,false,false,true);
        fprintf('OK.\n');
    catch
        fprintf('Not OK.\n');
    end
end

%Check if it is possible to import an SBML model using libSBML
fprintf(' > Checking import of model in SBML format...  ');
if isOctave
    fprintf('Skipped, MATLAB only.\n')
else
    try
        evalc('importModel(xmlFile)');
        try
            libSBMLver=OutputSBML; % Only works in libSBML 5.17.0+
            fprintf('OK.\n');
        catch
            fprintf(['Not OK.\n\n'...
                'An older libSBML version was found, update to version 5.17.0 or higher\n'...
                'for a significant improvement of model import\n\n']);
        end
    catch
        fprintf(['Not OK\nTo import SBML models, download libSBML from\n'...
            'http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n']);
    end
end

%Check if it is possible to import an YAML model
fprintf(' > Checking import of model in YAML format...  ');
try
    readYaml(ymlFile,true);
    fprintf('OK.\n');
catch
    fprintf('Not OK.\n');
end

fprintf('\n=== Model solvers ===\n');
%Define values for keepSolver and workingSolvers, needed for solver
%functionality check
keepSolver=false;
workingSolvers='';
%Get current solver. Set it to 'none', if it is not set
fprintf(' > Checking for solver in preferences...  ');
if ~ispref('RAVEN','solver')
    fprintf('None found\n')
    setRavenSolver('none');
    curSolv=getpref('RAVEN','solver');
else
    curSolv=getpref('RAVEN','solver');
    curSolv=regexprep(curSolv,'(_octave)|(_matlab)$','');
    fprintf('%s.\n',curSolv);
end

%Check if it is possible to solve an LP problem using different solvers
solver={'gurobi','glpk','cobra'};

for i=1:numel(solver)
    fprintf(' > Checking to solve an LP problem with %s...  ',solver{i});
    load(matFile);
    try
        setRavenSolver(solver{i});
        solveLP(emptyModel);
        workingSolvers=strcat(workingSolvers,';',solver{i});
        fprintf('OK.\n');
        if strcmp(curSolv,solver{i})
            keepSolver=true;
        end
    catch
        fprintf('Not OK.\n');
    end
end

fprintf(' > Solver set to...  ');
if keepSolver
    %The solver set in curSolv is functional, so the settings are restored
    %to the ones which were set before running checkInstallation
    setRavenSolver(curSolv);
    fprintf('%s.\n',curSolv);
elseif ~isempty(workingSolvers)
    %There are working solvers, but the none of them is the solver defined
    %by curSolv. The first working solver is therefore set as RAVEN solver
    workingSolvers=regexprep(workingSolvers,'^;','');
    workingSolvers=regexprep(workingSolvers,';.+$','');
    %Only one working solver should be left by now in workingSolvers
    setRavenSolver(workingSolvers);
    fprintf('%s.\n',workingSolvers);
else
    %No functional solvers were found, so the setting is restored back to
    %original
    setRavenSolver(curSolv);
    fprintf(['WARNING: No working solver was found!\n'...
        'Install the solver, set it using setRavenSolver(''solverName'') and run checkInstallation again\n'...
        'Available solverName options are ''gurobi'' and ''cobra''\n\n']);
end

fprintf('Checking essential binary executables:\n');
fprintf('NOTE: Broken binary executables <strong>must be fixed</strong> before running RAVEN\n');

fprintf('\tBLAST+... ');
res=runtests('blastPlusTests.m','OutputDetail',0);
interpretResults(res);
fprintf('\tDIAMOND... ');
res=runtests('diamondTests.m','OutputDetail',0);
interpretResults(res);
fprintf('\tHMMER... ');
res=runtests('hmmerTests.m','OutputDetail',0);
interpretResults(res);
fprintf('Checking non-essential/development binary executables:\n');
fprintf('NOTE: Only fix these binaries if planning to use KEGG FTP dump files in getKEGGModelForOrganism\n');
fprintf('\tCD-HIT... ');
res=runtests('cdhitTests.m','OutputDetail',0);
interpretResults(res);
fprintf('\tMAFFT... ');
res=runtests('mafftTests.m','OutputDetail',0);
interpretResults(res);

fprintf('\n=== Binaries ===\n');
fprintf(' > Checking essential binary executables for model reconstruction:\n');

tryBinary('makeblastdb',...
        ['"' fullfile(ravenDir,'software','blast+',['makeblastdb' binEnd]) '" -version'],...
        'Not OK! getBlast() will not work.')

tryBinary('blastp',...
        ['"' fullfile(ravenDir,'software','blast+',['blastp' binEnd]) '" -version'],...
        'Not OK! getBlast() will not work.')

tryBinary('diamond',...
        ['"' fullfile(ravenDir,'software','diamond',['diamond' binEnd]) '" version'],...
        'Not OK! getDiamond() will not work.')

tryBinary('hmmsearch',...
        ['"' fullfile(ravenDir,'software','hmmer',['hmmsearch' binEnd]) '" -h'],...
        'Not OK! getKEGGModelForOrganism() will not work.') % It can still make
        % model by providing kegg ID, but message is simplified.

fprintf('\n > Checking non-essential/development binary executables:\n');
fprintf(' > Note: These are only required when using KEGG FTP dump files in getKEGGModelForOrganism.\n');

tryBinary('cd-hit',...
        ['"' fullfile(ravenDir,'software','cd-hit',['cd-hit' binEnd]) '" -h'],...
        'Not OK! getKEGGModelForOrganism() will not work.')

if ismac
    mafftPath = ['"' fullfile(ravenDir,'software','mafft','mafft-mac','mafft.bat') '" --version '];
elseif isunix
    mafftPath = ['"' fullfile(ravenDir,'software','mafft','mafft-linux64','mafft.bat') '" --version '];
elseif ispc
    mafftPath = ['"' fullfile(ravenDir,'software','mafft','mafft-win','mafft.bat') '" --version '];
end
tryBinary('mafft', mafftPath,'Not OK! getKEGGModelForOrganism() will not work.')

tryBinary('hmmbuild',...
        ['"' fullfile(ravenDir,'software','hmmer',['hmmbuild' binEnd]) '" -h'],...
        'Not OK! getKEGGModelForOrganism() will not work.')

fprintf('\n=== Compatibility ===\n');
fprintf('Checking whether RAVEN functions are non-redundant across %s path...  ',mORo);
checkFunctionUniqueness();

fprintf('\n*** checkInstallation complete ***\n\n');
end

function interpretResults(results)
if results.Failed==0 && results.Incomplete==0
    fprintf('OK\n');
else
    fprintf('Not OK! Download/compile the binary and rerun checkInstallation\n');
end
end

function tryBinary(binName,fullCommand,notOKmsg)
fprintf(' > > %s...  ',binName)
if isOctave && strcmp(binName,'mafft'); fprintf('\n'); end
try
    [res,tmp]=system(fullCommand);
    if strcmp(binName,'cd-hit') && res==1; res=0; end
    if res==0
        fprintf('OK.\n');
    else
        fprintf('%s\n',notOKmsg);
    end
catch
    fprintf('%s\n',notOKmsg);
end
end
end
