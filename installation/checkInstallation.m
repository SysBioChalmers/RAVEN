function checkInstallation()
% checkInstallation
%   The purpose of this function is to check if all necessary functions are
%   installed and working

%   Usage: checkInstallation()
%
%   Rasmus Agren, 2014-01-06
%

fprintf('*** RAVEN TOOLBOX v. 1.9\n');

%Check if RAVEN is in the path list
paths=textscan(path,'%s','delimiter', pathsep);
paths=paths{1};

%Get the RAVEN path
[ST I]=dbstack('-completenames');
[ravenDir,crap1,crap2]=fileparts(fileparts(ST(I).file));

if ismember(ravenDir,paths)
    fprintf('Checking if RAVEN is in the Matlab path... PASSED\n');
else
    fprintf('Checking if RAVEN is in the Matlab path... FAILED\n');
    addMe=input('\tWould you like to add the RAVEN directory to the path list? Y/N\n','s');
    if strcmpi(addMe,'y')
        addpath(ravenDir);
        savepath
    end
end

excelFile=fullfile(ravenDir,'tutorial','empty.xlsx');
xmlFile=fullfile(ravenDir,'tutorial','empty.xml');

%Check if it is possible to parse an Excel file
try
    importExcelModel(excelFile,false,false,true);
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... PASSED\n');
catch
    fprintf('Checking if it is possible to parse a model in Microsoft Excel format... FAILED\n');
    if ispc==false %Print info for UNIX/MacOS
        fprintf('\tThis functionality uses Microsoft Excel COM server, which works best for the Windows version of Matlab\n');
    end
end

%Check if it is possible to import an SBML model using libSBML
try
    smallModel=importModel(xmlFile);
    fprintf('Checking if it is possible to import an SBML model using libSBML... PASSED\n');
catch
    fprintf('Checking if it is possible to import an SBML model using libSBML... FAILED\n');
end

%Check if it is possible to solve a LP problem using Mosek
try
    setRavenSolver('mosek')
    solveLP(smallModel);
    fprintf('Checking if it is possible to solve a LP problem using Mosek... PASSED\n');
catch
    fprintf('Checking if it is possible to solve a LP problem using Mosek... FAILED\n');
end

try
    setRavenSolver('gurobi')
    solveLP(smallModel);
    fprintf('Checking if it is possible to solve a LP problem using Gurobi... PASSED\n');
catch
    fprintf('Checking if it is possible to solve a LP problem using Gurobi... FAILED\n');
end
end