function exportToExcelFormat(model,filename)
% exportToExcelFormat
%   Exports a model structure to the Microsoft Excel model format
%
%   model       a model structure
%   filename    file name of the Excel file. Both xls and xlsx are supported.
%               In order to preserve backward compatibility this could also
%               be only a path, in which case the model is exported to a set
%               of tab-delimited text files instead. See exportToTabDelimited
%               for details regarding that functionality
%
%   The resulting Excel file can be used with importExcelModel/SBMLFromExcel
%   for modelling or to generate a SBML file.
%
%   NOTE: No checks are made regarding the correctness of the model. Use
%         checkModelStruct to identify problems in the model structure
%
%   Usage: exportToExcelFormat(model,filename)
%
%   Rasmus Agren, 2014-01-07
%

[filePath, A, B]=fileparts(filename);

if ~any(filePath)
   filePath=pwd; 
end

%If a path was used call on exportToTabDelimited instead
if ~any(A) || ~any(B)
    exportToTabDelimited(model,filename);
    return;
end

%Adds the required classes to the Java path
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file)));
poiPATH=fullfile(ravenPath,'software','apache-poi');
javaaddpath(fullfile(poiPATH,'dom4j-1.6.1.jar'));
javaaddpath(fullfile(poiPATH,'poi-3.8-20120326.jar'));
javaaddpath(fullfile(poiPATH,'poi-ooxml-3.8-20120326.jar'));
javaaddpath(fullfile(poiPATH,'poi-ooxml-schemas-3.8-20120326.jar'));
javaaddpath(fullfile(poiPATH,'xmlbeans-2.3.0.jar'));

%Import required classes from Apache POI
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.hssf.usermodel.*;
import org.apache.poi.xssf.usermodel.*;
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.util.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

%If the folder doesn't exist then create it
if ~exist(filePath,'dir')
    mkdir(filePath);
end

%Remove the file if it already exist
if exist(filename,'file')
    delete(filename);
end

%Check that the file endings are correct
if ~strcmpi(B,'.XLS') && ~strcmpi(B,'.XLSX')
   dispEM('The file name must end in .xls or .xlsx'); 
end

%Construct equations
model.equations=constructEquations(model,model.rxns,true);

%Check if it should print genes
if isfield(model,'grRules');    
    %Also do some parsing here
    rules=model.grRules;
    rules=strrep(rules,'(','');
    rules=strrep(rules,')','');
    rules=strrep(rules,' and ',':');
    rules=strrep(rules,' or ',';');
else
    rules=[];
end

%Check if the model has default upper/lower bounds. This determines if
%those values should be printed or not
hasDefaultLB=false;
hasDefaultUB=false;
if isfield(model,'annotation')
    if isfield(model.annotation,'defaultLB')
       hasDefaultLB=true; 
    end
    if isfield(model.annotation,'defaultUB')
       hasDefaultUB=true; 
    end
end

%Create the workbook
if strcmpi(B,'.XLS')
    wb=HSSFWorkbook();
else
    wb=XSSFWorkbook();
end

%Add the RXNS sheet
s=wb.createSheet();
wb.setSheetName(0, 'RXNS');

%Create the header row
headers={'#';'ID';'NAME';'EQUATION';'EC-NUMBER';'GENE ASSOCIATION';'LOWER BOUND';'UPPER BOUND';'OBJECTIVE';'COMPARTMENT';'MIRIAM';'SUBSYSTEM';'REPLACEMENT ID'};
r=s.createRow(0);
for i=0:numel(headers)-1
    c=r.createCell(i);
    c.setCellValue(headers{i+1});    
end

%Then fill in the sheet
for i=1:numel(model.rxns)
    r=s.createRow(i);
    
    c=r.createCell(1);
    c.setCellValue(model.rxns{i});
    
    if isfield(model,'rxnNames')
        c=r.createCell(2);
        c.setCellValue(model.rxnNames{i});
    end

    c=r.createCell(3);
    c.setCellValue(model.equations{i});
    
    if isfield(model,'eccodes')
        c=r.createCell(4);
        c.setCellValue(model.eccodes{i});
    end
    
    if ~isempty(rules)
        c=r.createCell(5);
        c.setCellValue(rules{i});
    end
    
    if isfield(model,'lb')
        if hasDefaultLB==true
            if model.rev(i)==1
                %If reversible, print only if different than defaultLB
                if model.lb(i) ~= model.annotation.defaultLB
                    c=r.createCell(6);
                    c.setCellValue(model.lb(i));
                end
            else
                %If irreversible, print only for non-zero values
            	if model.lb(i)~=0
                    c=r.createCell(6);
                    c.setCellValue(model.lb(i));
            	end
            end
        else
            c=r.createCell(6);
            c.setCellValue(model.lb(i));
        end
    end
    if isfield(model,'ub')
        if hasDefaultUB==true
            if model.ub(i) ~= model.annotation.defaultUB
                c=r.createCell(7);
                c.setCellValue(model.ub(i));
            end
        else
            c=r.createCell(7);
            c.setCellValue(model.ub(i));
        end
    end
    if isfield(model,'c')
        if model.c(i)~=0
            c=r.createCell(8);
            c.setCellValue(model.c(i));
        end
    end
    if isfield(model,'rxnComps')
        c=r.createCell(9);
        c.setCellValue(model.comps{model.rxnComps(i)});
    end
    if isfield(model,'rxnMiriams')
       if ~isempty(model.rxnMiriams{i})
           toPrint=[];
           for j=1:numel(model.rxnMiriams{i}.name)
               toPrint=[toPrint strtrim(model.rxnMiriams{i}.name{j}) ':' strtrim(model.rxnMiriams{i}.value{j}) ';']; 
           end
           c=r.createCell(10);
           c.setCellValue(toPrint(1:end-1));
       end 
    end
    if isfield(model,'subSystems')
        c=r.createCell(11);
        c.setCellValue(model.subSystems{i});
    end
end

%Add the METS sheet
s=wb.createSheet();
wb.setSheetName(1, 'METS');

%Create the header row
headers={'#';'ID';'NAME';'UNCONSTRAINED';'MIRIAM';'COMPOSITION';'InChI';'COMPARTMENT';'REPLACEMENT ID'};
r=s.createRow(0);
for i=0:numel(headers)-1
    c=r.createCell(i);
    c.setCellValue(headers{i+1});    
end

for i=1:numel(model.mets)
    r=s.createRow(i);
    
    c=r.createCell(1);
    c.setCellValue([model.metNames{i} '[' model.comps{model.metComps(i)} ']']);
    
    if isfield(model,'metNames')
        c=r.createCell(2);
        c.setCellValue(model.metNames{i});
    end
    
    if isfield(model,'unconstrained')
        if model.unconstrained(i)~=0
            c=r.createCell(3);
            c.setCellValue(true);
        end
    end
    
    if isfield(model,'metMiriams')
       if ~isempty(model.metMiriams{i})
           toPrint=[];
           for j=1:numel(model.metMiriams{i}.name)
               toPrint=[toPrint strtrim(model.metMiriams{i}.name{j}) ':' strtrim(model.metMiriams{i}.value{j}) ';']; 
           end
           c=r.createCell(4);
           c.setCellValue(toPrint(1:end-1));
       end 
    end
    
    if isfield(model,'metFormulas')
        c=r.createCell(5);
        c.setCellValue(model.metFormulas{i});
    end
    
    if isfield(model,'inchis')
        c=r.createCell(6);
        c.setCellValue(model.inchis{i});
    end
    
    if isfield(model,'metComps')
        c=r.createCell(7);
        c.setCellValue(model.comps{model.metComps(i)});
    end
    
    c=r.createCell(8);
    c.setCellValue(model.mets{i});
end

%Add the COMPS sheet
s=wb.createSheet();
wb.setSheetName(2, 'COMPS');

%Create the header row
headers={'#';'ABBREVIATION';'NAME';'INSIDE';'MIRIAM'};
r=s.createRow(0);
for i=0:numel(headers)-1
    c=r.createCell(i);
    c.setCellValue(headers{i+1});    
end

for i=1:numel(model.comps)
    r=s.createRow(i);
    
    c=r.createCell(1);
    c.setCellValue(model.comps{i});
    
    if isfield(model,'compNames')
        c=r.createCell(2);
        c.setCellValue(model.compNames{i});
    end
    
    if isfield(model,'compOutside')
        c=r.createCell(3);
        c.setCellValue(model.compOutside{i});
    end
    
    if isfield(model,'compMiriams')
       if ~isempty(model.compMiriams{i})
           toPrint=[];
           for j=1:numel(model.compMiriams{i}.name)
               toPrint=[toPrint strtrim(model.compMiriams{i}.name{j}) ':' strtrim(model.compMiriams{i}.value{j}) ';']; 
           end
           c=r.createCell(4);
           c.setCellValue(toPrint(1:end-1));
       end 
    end
end

%Add the GENES sheet
addedGeneSheet=-1; %This is to get the sheet numbering right if no genes are added
if isfield(model,'genes')
    s=wb.createSheet();
    wb.setSheetName(3, 'GENES');
    addedGeneSheet=0;

    %Create the header row
    headers={'#';'NAME';'MIRIAM';'SHORT NAME';'COMPARTMENT'};
    r=s.createRow(0);
    for i=0:numel(headers)-1
        c=r.createCell(i);
        c.setCellValue(headers{i+1});    
    end

    for i=1:numel(model.genes)
        r=s.createRow(i);

        c=r.createCell(1);
        c.setCellValue(model.genes{i});
        
       if isfield(model,'geneMiriams')
           if ~isempty(model.geneMiriams{i})
               toPrint=[];
               for j=1:numel(model.geneMiriams{i}.name)
                   toPrint=[toPrint strtrim(model.geneMiriams{i}.name{j}) ':' strtrim(model.geneMiriams{i}.value{j}) ';']; 
               end
               c=r.createCell(2);
               c.setCellValue(toPrint(1:end-1));
           end 
       end
       
       if isfield(model,'geneComps')
           c=r.createCell(4);
           c.setCellValue(model.comps{model.geneComps(i)});
       end
    end
end

%Add the MODEL sheet
s=wb.createSheet();
wb.setSheetName(4-addedGeneSheet, 'MODEL');

%Create the header row
headers={'#';'ID';'DESCRIPTION';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'TAXONOMY';'NOTES'};
r=s.createRow(0);
for i=0:numel(headers)-1
    c=r.createCell(i);
    c.setCellValue(headers{i+1});    
end

%Add some default stuff if needed
if ~isfield(model,'annotation')
   model.annotation.familyName='Agren';
   model.annotation.givenName='Rasmus';
   model.annotation.email='rasmus.j.agren@gmail.com';
   model.annotation.organization='Chalmers University of Technology';
end

r=s.createRow(1);
a=model.annotation;

if isfield(model,'id')
    c=r.createCell(1);
    c.setCellValue(model.id);
end
if isfield(model,'description')
    c=r.createCell(2);
    c.setCellValue(model.description);
end
if isfield(a,'defaultLB')
    c=r.createCell(3);
    c.setCellValue(a.defaultLB);
end
if isfield(a,'defaultUB')
    c=r.createCell(4);
    c.setCellValue(a.defaultUB);
end
if isfield(a,'givenName')
    c=r.createCell(5);
    c.setCellValue(a.givenName);
end
if isfield(a,'familyName')
    c=r.createCell(6);
    c.setCellValue(a.familyName);
end
if isfield(a,'email')
    c=r.createCell(7);
    c.setCellValue(a.email);
end
if isfield(a,'organization')
    c=r.createCell(8);
    c.setCellValue(a.organization);
end
if isfield(a,'taxonomy')
    c=r.createCell(9);
    c.setCellValue(a.taxonomy);
end
if isfield(a,'note')
    c=r.createCell(10);
    c.setCellValue(a.note);
end
        
%Open the output stream
out = FileOutputStream(filename);
wb.write(out);
out.close();
end
