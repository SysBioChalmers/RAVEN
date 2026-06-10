function SBMLFromExcel(fileName, outputFileName,toCOBRA,printWarnings)
% SBMLFromExcel  Convert a model in the Excel format to SBML.
%
% For a detailed description of the file format, see the supplied manual.
%
% Parameters
% ----------
% fileName : char
%     the Excel file.
% outputFileName : char
%     the SBML file.
% toCOBRA : logical, optional
%     true if the model should be saved in COBRA Toolbox format. Only
%     limited support at the moment (default false).
% printWarnings : logical, optional
%     true if warnings about model issues should be reported (default
%     true).
%
% Examples
% --------
%     SBMLFromExcel(fileName, outputFileName, toCOBRA, printWarnings);
%
% Notes
% -----
% This is just a wrapper function for importExcelModel, printModelStats and
% exportModel. Use those functions directly for greater control.
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
