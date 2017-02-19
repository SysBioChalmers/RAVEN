% loadWorkbook
%   Loads an Excel file into a Workbook object using the Java library Apache POI
%
%   fileName    name of the Excel file. If it doesn't exist it will be
%               created
%   createEmpty true if an empty workbook should be created if the file
%               didn't exist (opt, default false)
%
%   workbook    Workbook object representing the Excel file
%
%   Usage: workbook=loadWorkbook(fileName,createEmpty)
%
%   Rasmus Agren, 2015-08-19
%
function workbook=loadWorkbook(fileName,createEmpty)
    if nargin<2
        createEmpty=false;
    end
    
    %Adds the required classes to the static Java path
    addJavaPaths();

    %Import required classes from Apache POI
    import org.apache.poi.ss.usermodel.*;
    import org.apache.poi.ss.util.*;
    import java.io.FileInputStream;
    import org.apache.poi.hssf.usermodel.*;
    import org.apache.poi.xssf.usermodel.*;
    
    %Check if the file exists
    if ~exist(fileName,'file')
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