% loadSheet
%   Loads an Excel sheet into a cell matrix using the Java library Apache POI
%
%   workbook    Workbook object representing the Excel file
%   sheet       name of the sheet (opt, default first sheet)
%
%   raw         cell array with the data in the sheet
%   flag        0 if everything worked, -1 if it didn't
%
%   Usage: [raw, flag]=loadSheet(workbook, sheet)

function [raw, flag]=loadSheet(workbook, sheet)
if nargin<2
    sheet=[];
end
flag=0;
raw={};

if any(sheet)
    sh=workbook.getSheet(sheet);
else
    sh=workbook.getSheetAt(0);
end
if isempty(sh)
    flag=-1;
    return;
end

lastRow=sh.getLastRowNum();
wasEmpty=false(lastRow+1,1);
raw=cell(lastRow+1,0); %Allocate space for the cell array. The number of columns isn't know yet, as it's saved row-by-row
for i=0:lastRow
    row=sh.getRow(i);
    %Sometimes the last rows only contain formatting (or some other weird
    %Excel thing). Ignore such empty rows. Note the +1 to deal with that
    %Matlab indexing starts at 1
    if isempty(row)
        wasEmpty(i+1)=true;
        continue;
    end
    lastCol=row.getLastCellNum();
    
    %Adjust the size of the cell array if needed
    if (lastCol+1)>size(raw,2)
        raw=[raw cell(lastRow+1,lastCol+1-size(raw,2))];
    end
    
    %Loop over the columns
    for j=0:lastCol
        c=row.getCell(j,row.RETURN_BLANK_AS_NULL);
        if ~isempty(c)
            %Then decide how to save it depending on the type. First check
            %if it's a formula. If so the cached value should be used
            if c.getCellType()==c.CELL_TYPE_FORMULA
                type=c.getCachedFormulaResultType();
            else
                type=c.getCellType();
            end
            
            switch type
                case c.CELL_TYPE_STRING
                    raw{i+1,j+1}=char(c.getRichStringCellValue().getString());
                case c.CELL_TYPE_NUMERIC
                    raw{i+1,j+1}=c.getNumericCellValue();
                case c.CELL_TYPE_BOOLEAN
                    raw{i+1,j+1}=c.getBooleanCellValue();
            end
        end
    end
end

%Remove empty rows
raw(wasEmpty,:)=[];
end
