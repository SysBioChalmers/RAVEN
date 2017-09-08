function checkInstallation()
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working. It also checks whether there ary any functions
%   with overlapping names between RAVEN and other toolboxes or
%   user-defined functions, which are accessible from Matlab pathlist
%
%   Usage: checkInstallation()
%
%	Simonas Marcisauskas, 2017-09-08
%

fprintf('*** RAVEN TOOLBOX v. 1.9\n');

keepSolver=false;
lastWorking='';

%Check if RAVEN is in the path list
paths=textscan(path,'%s','delimiter', pathsep);
paths=paths{1};

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

% get current solver
if ~ispref('RAVEN','solver')
	fprintf('Solver found in preferences... NONE\n');
else
	curSolv=getpref('RAVEN','solver');
	fprintf(['Solver found in preferences... ',curSolv,'\n']);
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

%Check if it is possible to parse an Excel file
try
    importExcelModel(excelFile,false,false,true);
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... PASSED\n');
catch
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... FAILED\n');
end

%Check if it is possible to import an SBML model using libSBML
try
    smallModel=importModel(xmlFile);
    fprintf('Checking if it is possible to import an SBML model using libSBML... PASSED\n');
catch
    fprintf('Checking if it is possible to import an SBML model using libSBML... FAILED\nTo import SBML models, download libSBML from http://sbml.org/Software/libSBML/Downloading_libSBML and add to MATLAB path\n');
end

%Check if it is possible to solve an LP problem using different solvers
solver={'mosek','gurobi'};

for i=1:numel(solver)
    try
        setRavenSolver(solver{i});
        solveLP(smallModel);
        lastWorking=solver{i};
        fprintf(['Checking if it is possible to solve an LP problem using ',solver{i},'... PASSED\n']);
        if and(exist('curSolv','var'),strcmp(curSolv,solver{i}))
            keepSolver=true;
        end
    catch
        fprintf(['Checking if it is possible to solve an LP problem using ',solver{i},'... FAILED\n']);
    end
end

if keepSolver
    setRavenSolver(curSolv);
elseif ~isempty(lastWorking)
    setRavenSolver(lastWorking);
end

if ~exist('curSolv','var')
	fprintf(['Preferred solver... NEW\nSolver saved as preference... ',lastWorking,'\n']);
elseif keepSolver
	fprintf(['Preferred solver... KEPT\nSolver saved as preference... ',curSolv,'\n']);
else
	fprintf(['Preferred solver... CHANGED\nSolver saved as preference... ',lastWorking,'\n']);
end

fprintf(['Checking the uniqueness of RAVEN functions across Matlab path...\n']);
checkFunctionUniqueness();

end
