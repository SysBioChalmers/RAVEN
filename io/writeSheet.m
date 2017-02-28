function wb=writeSheet(wb,sheetName,sheetPosition,captions,units,raw)
% writeSheet
%   Writes a cell matrix to an Excel sheet into using the Java library Apache POI
%
%   workbook        Workbook object representing the Excel file
%   sheetName       name of the sheet
%   sheetPosition   0-based position of the sheet
%   captions        cell array of captions (opt)
%   units            WRITE INFO
%   raw             cell array with the data in the sheet
%
%   Usage: wb=writeSheet(wb,sheetName,sheetPosition,captions,units,raw)
%
%   Rasmus Agren, 2017-02-27
%

%Adds the required classes to the static Java path if not already added
addJavaPaths();

%Import required classes from Apache POI.
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.util.*;
import org.apache.poi.hssf.usermodel.*;
import org.apache.poi.xssf.usermodel.*;
    
%Set default style object and bold style object
defaultStyle = wb.createCellStyle();
idx = wb.getCreationHelper().createDataFormat().getFormat('##0.00');
defaultStyle.setDataFormat(idx);

boldFont=wb.createFont();
boldFont.setBoldweight(boldFont.BOLDWEIGHT_BOLD);
boldStyle=defaultStyle.clone();
boldStyle.setFont(boldFont);

s=wb.createSheet();
wb.setSheetName(sheetPosition, sheetName);

%Create the header row and units
counter=0;
if ~isempty(captions)
    r=s.createRow(counter);
    for i=0:numel(captions)-1
        c=r.createCell(i);
        c.setCellValue(captions{i+1});
        c.setCellStyle(boldStyle);
    end
    counter=counter+1;
end
if ~isempty(units)
    r=s.createRow(counter);
    for i=0:numel(units)-1
        c=r.createCell(i);
        content=units{i+1};
        if iscell(content) && numel(content)==2
            c.setCellValue(content{1});

            %This means that the cell should have a list of allowed
            %values
            dvHelper=XSSFDataValidationHelper(s);

            cellRange=CellRangeAddressList(1,1,i,i); %Units are on the second line
            dvConstraint=dvHelper.createExplicitListConstraint(toValid(content{2}));
            validation=dvHelper.createValidation(dvConstraint, cellRange);
            validation.setShowErrorBox(true);
            s.addValidationData(validation);
        else
            c.setCellValue(content{1});
        end
        c.setCellStyle(boldStyle);
    end
    counter=counter+1;
end

%Loop through and fill in the values
for i=0:size(raw,1)-1
    r=s.createRow(i+counter);
    
    for j=0:size(raw,2)-1
        if ~isempty(raw{i+1,j+1})
            c=r.createCell(j);
            content=raw{i+1,j+1};
            if iscell(content) && numel(content)==2
                c.setCellValue(content{1});
                
                %This means that the cell should have a list of allowed
                %values
                dvHelper=XSSFDataValidationHelper(s);
                cellRange=CellRangeAddressList(i+1,i+1,j,j);
                dvConstraint=dvHelper.createExplicitListConstraint(toValid(content{2}));
                validation=dvHelper.createValidation(dvConstraint, cellRange);
                validation.setShowErrorBox(true);
                s.addValidationData(validation);
            else
                c.setCellValue(content);
            end
            c.setCellStyle(defaultStyle);
        end
    end
end

%Resize columns
for i=0:size(raw,2)-1
    s.autoSizeColumn(i);
end

%Add freeze panes
if counter>0
    s.createFreezePane(0,counter,0,counter);
end
end

function I=toValid(J)
    I=cell(numel(J),1);
    for i=1:numel(J)
        if ischar(J{i})
            I(i)=J(i);
        else
            if islogical(J{i})
               if J{i}==true
                    I{i}='TRUE';
               else
                    I{i}='FALSE';
               end
            else
                %Other, most likely numbers
                I{i}=num2str(J{i});
            end
        end
    end
end