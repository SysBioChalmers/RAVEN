function workbook=loadWorkbook(fileName,varargin)
% loadWorkbook  Load an Excel file into a Workbook object.
%
% Loads an Excel file into a Workbook object using the Java library Apache
% POI.
%
% Parameters
% ----------
% fileName : char
%     name of the Excel file. If it doesn't exist it will be created.
%
% Name-Value Arguments
% --------------------
% createEmpty : logical
%     true if an empty workbook should be created if the file didn't exist
%     (default false).
%
% Returns
% -------
% workbook : Workbook
%     Workbook object representing the Excel file.
%
% Examples
% --------
%     workbook = loadWorkbook(fileName, createEmpty);

p=parseRAVENargs(varargin, {'createEmpty',false});
createEmpty=p.createEmpty;

%Check if the user has MATLAB Text Analytics Toolbox installed, as it comes
%with its own conflicting version of the required Apache POI files
if exist('vaderSentimentScores.m')== 2
    error(['MATLAB Text Analytics Toolbox found. This should be uninstalled ' ...
           'if you want to read/write Excel files. See RAVEN GitHub Issues '...
           'page for instructions.'])    
end
    
%Adds the required classes to the static Java path if not already added
addJavaPaths();

%Import required classes from Apache POI
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.util.*;
import java.io.FileInputStream;
import org.apache.poi.hssf.usermodel.*;
import org.apache.poi.xssf.usermodel.*;

%Check if the file exists
if ~isfile(fileName)
    if createEmpty==false
        EM='The Excel file could not be found';
        dispEM(EM);
    else
        %Create an empty workbook
        [~,~,I]=fileparts(fileName);
        if strcmpi(I,'.xls')
            workbook=HSSFWorkbook();
        else
            if strcmpi(I,'.xlsx')
                workbook=XSSFWorkbook();
            else
                EM='The file name must end in .xls or .xlsx';
                dispEM(EM);
            end
        end
    end
else
    %Opens the workbook. The input stream is needed since it will otherwise
    %keep the file open
    is=FileInputStream(getFullPath(fileName));
    workbook=WorkbookFactory.create(is);
    is.close();
end
end
