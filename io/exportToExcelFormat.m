function exportToExcelFormat(model,filename)
% exportToExcelFormat
%   Exports a model structure to the Microsoft Excel model format
%
%   model       a model structure
%   filename    file name of the Excel file. Only xlsx format is supported.
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
%   Simonas Marcisauskas, 2017-08-25
%

[~, A, B]=fileparts(filename);

%If a path was used call on exportToTabDelimited instead
if ~any(A) || ~any(B)
    exportToTabDelimited(model,filename);
    return;
end

if ~strcmpi(B,'.xlsx')
    EM='As of RAVEN version 1.9, only export to xlsx format is supported';
    dispEM(EM);
end

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

%Remove the output file if it already exists
if exist(filename,'file')
    delete(filename);
end

%Load an empty workbook
wb=loadWorkbook(filename,true);

%Construct equations
model.equations=constructEquations(model,model.rxns,true);

%Check if it should print genes
if isfield(model,'grRules')
    rules=model.grRules;
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

%Add the RXNS sheet

%Create the header row
headers={'#';'ID';'NAME';'EQUATION';'EC-NUMBER';'GENE ASSOCIATION';'LOWER BOUND';'UPPER BOUND';'OBJECTIVE';'COMPARTMENT';'MIRIAM';'SUBSYSTEM';'REPLACEMENT ID';'NOTE';'REFERENCE';'CONFIDENCE SCORE'};

%Add empty comments
emptyColumn=cell(numel(model.rxns),1);
rxnSheet=emptyColumn;

%Add the model fields
rxnSheet=[rxnSheet model.rxns];

if isfield(model,'rxnNames')
    rxnSheet=[rxnSheet model.rxnNames];
else
    rxnSheet=[rxnSheet emptyColumn];
end

rxnSheet=[rxnSheet model.equations];

if isfield(model,'eccodes')
    rxnSheet=[rxnSheet model.eccodes];
else
    rxnSheet=[rxnSheet emptyColumn];
end

if ~isempty(rules)
    rxnSheet=[rxnSheet rules];
else
    rxnSheet=[rxnSheet emptyColumn];
end

lb=emptyColumn;
ub=emptyColumn;
objective=emptyColumn;
rxnMiriams=emptyColumn;

for i=1:numel(model.rxns)
    if isfield(model,'lb')
        if hasDefaultLB==true
            if model.rev(i)==1
                %If reversible, print only if different than defaultLB
                if model.lb(i) ~= model.annotation.defaultLB
                    lb{i}=model.lb(i);
                end
            else
                %If irreversible, print only for non-zero values
                if model.lb(i)~=0
                    lb{i}=model.lb(i);
                end
            end
        else
            lb{i}=model.lb(i);
        end
    end
    
    if isfield(model,'ub')
        if hasDefaultUB==true
            if model.ub(i) ~= model.annotation.defaultUB
                ub{i}=model.ub(i);
            end
        else
            ub{i}=model.ub(i);
        end
    end
    
    if isfield(model,'c')
        if model.c(i)~=0
            objective{i}=model.c(i);
        end
    end
    
    if isfield(model,'rxnMiriams')
       if ~isempty(model.rxnMiriams{i})
           toPrint=[];
           for j=1:numel(model.rxnMiriams{i}.name)
               toPrint=[toPrint strtrim(model.rxnMiriams{i}.name{j}) '/' strtrim(model.rxnMiriams{i}.value{j}) ';'];
           end
           rxnMiriams{i}=toPrint(1:end-1);
       end
    end
end

rxnSheet=[rxnSheet lb];
rxnSheet=[rxnSheet ub];
rxnSheet=[rxnSheet objective];

if isfield(model,'rxnComps')
    rxnSheet=[rxnSheet model.comps(model.rxnComps)];
else
    rxnSheet=[rxnSheet emptyColumn];
end

rxnSheet=[rxnSheet rxnMiriams];

if isfield(model,'subSystems')
    rxnSheet=[rxnSheet model.subSystems];
else
    rxnSheet=[rxnSheet emptyColumn];
end

%For REPLACEMENT ID which isn't in the model
rxnSheet=[rxnSheet emptyColumn];

if isfield(model,'rxnNotes')
    rxnSheet=[rxnSheet model.rxnNotes];
else
    rxnSheet=[rxnSheet emptyColumn];
end

if isfield(model,'rxnReferences')
    rxnSheet=[rxnSheet model.rxnReferences];
else
    rxnSheet=[rxnSheet emptyColumn];
end

if isfield(model,'rxnConfidenceScores')
    rxnSheet=[rxnSheet model.rxnConfidenceScores];
else
    rxnSheet=[rxnSheet emptyColumn];
end

wb=writeSheet(wb,'RXNS',0,headers,[],rxnSheet);

headers={'#';'ID';'NAME';'UNCONSTRAINED';'MIRIAM';'COMPOSITION';'InChI';'COMPARTMENT';'REPLACEMENT ID';'CHARGE'};

metSheet=cell(numel(model.mets),numel(headers));

for i=1:numel(model.mets)
    metSheet{i,2}=[model.metNames{i} '[' model.comps{model.metComps(i)} ']'];

    if isfield(model,'metNames')
        metSheet(i,3)=model.metNames(i);
    end

    if isfield(model,'unconstrained')
        if model.unconstrained(i)~=0
            metSheet{i,4}=true;
        end
    end

    if isfield(model,'metMiriams')
       if ~isempty(model.metMiriams{i})
           toPrint=[];
           for j=1:numel(model.metMiriams{i}.name)
               toPrint=[toPrint strtrim(model.metMiriams{i}.name{j}) '/' strtrim(model.metMiriams{i}.value{j}) ';'];
           end
           metSheet{i,5}=toPrint(1:end-1);
       end
    end

    if isfield(model,'metFormulas')
        metSheet(i,6)=model.metFormulas(i);
    end

    if isfield(model,'inchis')
        metSheet(i,7)=model.inchis(i);
    end

    if isfield(model,'metComps')
        metSheet(i,8)=model.comps(model.metComps(i));
    end

    metSheet(i,9)=model.mets(i);

    if isfield(model,'metCharge')
        metSheet{i,10}=model.metCharge(i);
    end
end

wb=writeSheet(wb,'METS',1,headers,[],metSheet);

%Add the COMPS sheet

%Create the header row
headers={'#';'ABBREVIATION';'NAME';'INSIDE';'MIRIAM'};

compSheet=cell(numel(model.comps),numel(headers));

for i=1:numel(model.comps)
    compSheet(i,2)=model.comps(i);

    if isfield(model,'compNames')
        compSheet(i,3)=model.compNames(i);
    end

    if isfield(model,'compOutside')
        compSheet(i,4)=model.compOutside(i);
    end

    if isfield(model,'compMiriams')
       if ~isempty(model.compMiriams{i})
           toPrint=[];
           for j=1:numel(model.compMiriams{i}.name)
               toPrint=[toPrint strtrim(model.compMiriams{i}.name{j}) '/' strtrim(model.compMiriams{i}.value{j}) ';'];
           end
           compSheet{i,5}=toPrint(1:end-1);
       end
    end
end

wb=writeSheet(wb,'COMPS',2,headers,[],compSheet);

%Add the GENES sheet
if isfield(model,'genes')
    %Create the header row
    headers={'#';'NAME';'MIRIAM';'SHORT NAME';'COMPARTMENT'};
    
    geneSheet=cell(numel(model.genes),numel(headers));

    for i=1:numel(model.genes)
       geneSheet(i,2)=model.genes(i);

       if isfield(model,'geneMiriams')
           if ~isempty(model.geneMiriams{i})
               toPrint=[];
               for j=1:numel(model.geneMiriams{i}.name)
                   toPrint=[toPrint strtrim(model.geneMiriams{i}.name{j}) '/' strtrim(model.geneMiriams{i}.value{j}) ';'];
               end
               geneSheet{i,3}=toPrint(1:end-1);
           end
       end
       if isfield(model,'geneShortNames')
           geneSheet(i,4)=model.geneShortNames(i);
       end
       if isfield(model,'geneComps')
           geneSheet(i,5)=model.comps(model.geneComps(i));
       end
    end
    
    wb=writeSheet(wb,'GENES',3,headers,[],geneSheet);
end

%Add the MODEL sheet

%Create the header row
headers={'#';'ID';'DESCRIPTION';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'TAXONOMY';'NOTES'};

modelSheet=cell(1,numel(headers));

%Add some default stuff if needed
if ~isfield(model,'annotation')
   model.annotation.familyName='Agren';
   model.annotation.givenName='Rasmus';
   model.annotation.email='rasmus.agren@scilifelab.se';
   model.annotation.organization='Chalmers University of Technology';
end

if isfield(model,'id')
    modelSheet{1,2}=model.id;
end
if isfield(model,'description')
    modelSheet{1,3}=model.description;
end
if isfield(model.annotation,'defaultLB')
    modelSheet{1,4}=model.annotation.defaultLB;
end
if isfield(model.annotation,'defaultUB')
    modelSheet{1,5}=model.annotation.defaultUB;
end
if isfield(model.annotation,'givenName')
    modelSheet{1,6}=model.annotation.givenName;
end
if isfield(model.annotation,'familyName')
    modelSheet{1,7}=model.annotation.familyName;
end
if isfield(model.annotation,'email')
    modelSheet{1,8}=model.annotation.email;
end
if isfield(model.annotation,'organization')
    modelSheet{1,9}=model.annotation.organization;
end
if isfield(model.annotation,'taxonomy')
    modelSheet{1,10}=model.annotation.taxonomy;
end
if isfield(model.annotation,'note')
    modelSheet{1,11}=model.annotation.note;
end

if isfield(model,'genes')
    wb=writeSheet(wb,'MODEL',4,headers,[],modelSheet);
else
    wb=writeSheet(wb,'MODEL',3,headers,[],modelSheet);
end

%Open the output stream
out = FileOutputStream(filename);
wb.write(out);
out.close();
end
