function SBMLFromExcel(fileName, outputFileName,toCOBRA,printWarnings)
% SBMLFromExcel
%   Converts a model in the Excel format to SBML
%
%   fileName        the Excel file
%   outputFileName  the SBML file
%   toCOBRA         true if the model should be saved in COBRA Toolbox
%                   format. Only limited support at the moment (opt,
%                   default false)
%   printWarnings   true if warnings about model issues should be reported
%                   (opt, default true)
%
%   For a detailed description of the file format, see the supplied manual.
%
%   Usage: SBMLFromExcel(fileName,outputFileName,toCOBRA,printWarnings)
%
%   NOTE: This is just a wrapper function for importExcelModel, printModelStats
%   and exportModel. Use those functions directly for greater control.
fileName=char(fileName);
outputFileName=char(outputFileName);
if nargin<3
    toCOBRA=false;
end
if nargin<4
    printWarnings=true;
end

model=importExcelModel(fileName,false,printWarnings);
printModelStats(model,printWarnings,false);
exportModel(model,outputFileName,toCOBRA,true);
end
