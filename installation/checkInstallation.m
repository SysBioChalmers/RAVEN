function checkInstallation()
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there are any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from Matlab pathlist
%
%   Usage: checkInstallation()
%
%	Simonas Marcisauskas, 2019-10-04
%

%Check if RAVEN is in the Matlab path list
paths=textscan(path,'%s','delimiter', pathsep);
paths=paths{1};

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

%Print the RAVEN version if it is not the development version
if exist(fullfile(ravenDir,'version.txt'), 'file') == 2
    fprintf(['\n*** THE RAVEN TOOLBOX v.' fgetl(fopen(fullfile(ravenDir,'version.txt'))) ' ***\n\n']);
else
    fprintf('\n*** THE RAVEN TOOLBOX - DEVELOPMENT VERSION ***\n\n');
end

if ismember(ravenDir,paths)
    fprintf('Checking if RAVEN is on the Matlab path... PASSED\n');
else
    fprintf('Checking if RAVEN is on the Matlab path... FAILED\n');
    addMe=input('Would you like to add the RAVEN directory to the path list? Y/N\n','s');
    if strcmpi(addMe,'y')
        subpath=regexp(genpath(ravenDir),pathsep,'split'); % Lists all subdirectories
        pathsToKeep=cellfun(@(x) isempty(strfind(x,'.git')),subpath) & cellfun(@(x) isempty(strfind(x,'doc')),subpath);
        addpath(strjoin(subpath(pathsToKeep),pathsep));
        savepath
    end
end

%Adds the required classes to the static Java path if not already added
addJavaPaths();

excelFile=fullfile(ravenDir,'tutorial','empty.xlsx');
xmlFile=fullfile(ravenDir,'tutorial','empty.xml');
matFile=fullfile(ravenDir,'tutorial','empty.mat');

%Check if it is possible to parse an Excel file
try
    importExcelModel(excelFile,false,false,true);
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... PASSED\n');
catch
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... FAILED\n');
end

%Check if it is possible to import an SBML model using libSBML
try
    importModel(xmlFile);
    try
        libSBMLver=OutputSBML; % Only works in libSBML 5.17.0+
        fprintf('Checking if it is possible to import an SBML model using libSBML... PASSED\n');
    catch
        fprintf(['Checking if it is possible to import an SBML model using libSBML... PASSED\n\n'...
            'An older libSBML version was found, update to version 5.17.0 or higher\n'...
            'for a significant improvement of model import.\n\n']);
    end
catch
    fprintf('Checking if it is possible to import an SBML model using libSBML... FAILED\nTo import SBML models, download libSBML from http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n');
end

%Define values for keepSolver and workingSolvers, needed for solver
%functionality check
keepSolver=false;
workingSolvers='';
% Get current solver. Set it to 'none', if it is not set;
if ~ispref('RAVEN','solver')
    fprintf('Solver found in preferences... NONE\n');
    setRavenSolver('none');
    curSolv=getpref('RAVEN','solver');
else
    curSolv=getpref('RAVEN','solver');
    fprintf(['Solver found in preferences... ',curSolv,'\n']);
end

%Check if it is possible to solve an LP problem using different solvers
solver={'gurobi','mosek','cobra'};

for i=1:numel(solver)
    try
        setRavenSolver(solver{i});
        load(matFile);
        solveLP(empty);
        workingSolvers=strcat(workingSolvers,';',solver{i});
        fprintf(['Checking if it is possible to solve an LP problem using ',solver{i},'... PASSED\n']);
        if strcmp(curSolv,solver{i})
            keepSolver=true;
        end
    catch
        fprintf(['Checking if it is possible to solve an LP problem using ',solver{i},'... FAILED\n']);
    end
end

if keepSolver
    % The solver set in curSolv is functional, so the settings are restored
    % to the ones which were set before running checkInstallation;
    setRavenSolver(curSolv);
    fprintf(['Preferred solver... KEPT\nSolver saved as preference... ',curSolv,'\n\n']);
elseif ~isempty(workingSolvers)
    % There are working solvers, but the none of them is the solver defined
    % by curSolv. The first working solver is therefore set as RAVEN
    % solver;
    workingSolvers=regexprep(workingSolvers,'^;','');
    workingSolvers=regexprep(workingSolvers,';.+$','');
    % Only one working solver should be left by now in workingSolvers;
    setRavenSolver(workingSolvers);
    fprintf(['Preferred solver... NEW\nSolver saved as preference... ',workingSolvers,'\n\n']);
else
    % No functional solvers were found, so the setting is restored back to
    % original;
    setRavenSolver(curSolv);
    fprintf('WARNING: No working solver was found!\nInstall the solver, set it using setRavenSolver(''solverName'') and run checkInstallation again.\nAvailable solverName options are ''mosek'', ''gurobi'' and ''cobra''\n\n');
end

if ~ispc
    if ismac
        binEnd='.mac';
    elseif isunix
        binEnd='';
    end
    fprintf('Checking binary executables..\n');
    [res,~]=system(['"' fullfile(ravenDir,'software','blast+',['blastp' binEnd]) '"']);
    if res==1
        fprintf(['Checking blastp' binEnd '.. OK\n']);
    else
        fprintf(['Checking blastp' binEnd '.. Not OK! The binary must be recompiled from source before running RAVEN\n']);
    end
    [res,~]=system(['"' fullfile(ravenDir,'software','blast+',['makeblastdb' binEnd]) '"']);
    if res==1
        fprintf(['Checking makeblastdb' binEnd '.. OK\n']);
    else
        fprintf(['Checking makeblastdb' binEnd '.. Not OK! The binary must be recompiled from source before running RAVEN\n']);
    end
    [res,~]=system(['"' fullfile(ravenDir,'software','cd-hit',['cd-hit' binEnd]) '"']);
    if res==1
        fprintf(['Checking cd-hit' binEnd '.. OK\n']);
    else
        fprintf(['Checking cd-hit' binEnd '.. Not OK! The binary must be recompiled from source before running RAVEN\n']);
    end
    [res,~]=system(['"' fullfile(ravenDir,'software','hmmer',['hmmbuild' binEnd]) '"']);
    if res==1
        fprintf(['Checking hmmbuild' binEnd '.. OK\n']);
    else
        fprintf(['Checking hmmbuild' binEnd '.. Not OK! The binary must be recompiled from source before running RAVEN\n']);
    end
    [res,~]=system(['"' fullfile(ravenDir,'software','hmmer',['hmmsearch' binEnd]) '"']);
    if res==1
        fprintf(['Checking hmmsearch' binEnd '.. OK\n']);
    else
        fprintf(['Checking hmmsearch' binEnd '.. Not OK! The binary must be recompiled from source before running RAVEN\n']);
    end
    if ismac
        [res,~]=system(['"' fullfile(ravenDir,'software','mafft','mafft-mac','mafft.bat') '" --help ']);
    elseif unix
        [res,~]=system(['"' fullfile(ravenDir,'software','mafft','mafft-linux64','mafft.bat') '" --help ']);
    end
    if res==1
        fprintf('Checking mafft.bat.. OK\n\n');
    else
        fprintf('Checking mafft.bat.. Not OK! The binary must be recompiled from source before running RAVEN\n\n');
    end
end

fprintf('Checking the uniqueness of RAVEN functions across Matlab path...\n');
checkFunctionUniqueness();

end
