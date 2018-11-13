function blastStructure=getBlastFromExcel(models,blastFile,organismId)
% getBlastFromExcel
%   Retrieves gene homology information from Excel files. Used as
%   input to getModelFromHomology.
%
%   models          a cell array of model structures
%   blastFile       Excel file with homology information
%   organismId      the id of the organism of interest (as described in the
%                   Excel file)
%
%   blastStructure  structure containing the information in the Excel
%                   sheets.
%
%   The Excel file should contain a number of spreadsheets which in turn
%   contain the bidirectional homology measurements between the genes in the
%   organisms. The first and second column headers in each sheet is the
%   "to" and "from" model ids (as defined in models or for the new organism).
%   The entries should correspond to the gene names in those models. The third,
%   fourth, fifth, sixth and seventh columns represent the E-value, alignment
%   length, identity, bitscore and percentage of positive-scoring matches for
%   each measurement (captions should be "E-value", "Alignment length",
%   "Identity", "Bitscore" and "PPOS").
%
%   Usage: blastStructure=getBlastFromExcel(models,blastFile,organismId)
%
%   Eduard Kerkhoven, 2018-11-02
%

if ~(exist(blastFile,'file')==2)
    error('BLAST result file %s cannot be found',string(blastFile));
end

blastStructure=[];

%Get a list of model IDs
organisms=cell(numel(models)+1,1);
organisms{1}=organismId;
for i=1:numel(models)
    organisms{i+1}=models{i}.id;
end

%Get all the spreadsheets in the file
[type, sheets]=xlsfinfo(blastFile);

%Check if the file is a Microsoft Excel Spreadsheet
if ~any(regexp(type,'Excel Spreadsheet')
    EM='The file is not a Microsoft Excel Spreadsheet';
    dispEM(EM);
end

for i=1:numel(sheets)
    %Check if the sheet has the right header and deal with organisms that
    %are in "models"
    [values,dataSheet]=xlsread(blastFile,dataSheet{i});
    labels=dataSheet(1,:);
    if strcmpi(labels{3},'E-value') && strcmpi(labels{4},'Alignment length') ...
            && strcmpi(labels{5},'Identity') && strcmpi(labels{6},'Bitscore') ...
            && strcmpi(labels{7},'PPOS')
        %At least one of the organisms must have a model
        fromID=find(strcmpi(labels{1},organisms));
        toID=find(strcmpi(labels{2},organisms));
        %Check that the organism ids exist and that one of them is the
        %organism of interest
        if any(fromID) && any(toID) && (toID==1 || fromID==1)
            %Check that no gene ids are empty. This could for example be
            %the case if the gene names are wrongly formatted as numbers
            %instead of strings
            emptyNames=cellfun(@isempty,dataSheet(2:end,1)) | cellfun(@isempty,dataSheet(2:end,2));
            if any(emptyNames)
                if all(emptyNames)
                    EM=['Only empty gene names in sheet from ' organisms{fromID} ' to ' organisms{toID}];
                    dispEM(EM);
                else
                    EM=['Empty gene names in sheet from ' organisms{fromID} ' to ' organisms{toID} '. Ignoring genes with empty names'];
                    dispEM(EM,false);
                end
            end
            blastStructure(numel(blastStructure)+1).toId=organisms{toID};
            blastStructure(numel(blastStructure)).fromId=organisms{fromID};
            blastStructure(numel(blastStructure)).fromGenes=dataSheet(2:end,1);
            blastStructure(numel(blastStructure)).toGenes=dataSheet(2:end,2);
            blastStructure(numel(blastStructure)).evalue=values(:,1);
            blastStructure(numel(blastStructure)).aligLen=values(:,2);
            blastStructure(numel(blastStructure)).identity=values(:,3);
            blastStructure(numel(blastStructure)).bitscore=values(:,4);
            blastStructure(numel(blastStructure)).ppos=values(:,5);
            
            %Remove matches where any of the values is NaN. This would have
            %been done anyways in getModelFromHomology, but it's neater to
            %do it here
            I=isnan(blastStructure(end).evalue) | isnan(blastStructure(end).aligLen) | isnan(blastStructure(end).identity);
            blastStructure(end).fromGenes(I)=[];
            blastStructure(end).toGenes(I)=[];
            blastStructure(end).evalue(I)=[];
            blastStructure(end).aligLen(I)=[];
            blastStructure(end).identity(I)=[];
            blastStructure(end).bitscore(I)=[];
            blastStructure(end).ppos(I)=[];            
        else
            if isempty(toID) || isempty(fromID)
                EM=['The data in sheet ' sheets{i} ' has no corresponding model. Ignoring sheet'];
                dispEM(EM,false);
            else
                EM=['The data in sheet ' sheets{i} ' does not involve the organism of interest. Ignoring sheet'];
                dispEM(EM,false);
            end
        end
    else
        EM=['The data in sheet ' sheets{i} ' is not correctly formatted. Ignoring sheet'];
        dispEM(EM,false);
    end
end

end
