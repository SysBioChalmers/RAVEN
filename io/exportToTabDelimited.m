function exportToTabDelimited(model,path)
% exportToTabDelimited
%   Exports a model structure to a set of tab-delimited text files
%
%   model	a model structure
%   path	the path to export to. The resulting text files will be saved
%           under the names excelRxns.txt, excelMets.txt, excelGenes.txt,
%           excelModel.txt, and excelComps.txt
%
%   NOTE: This functionality was previously a part of exportToExcelFormat.
%         The naming of the resulting text files is to preserve backward
%         compatibility
%
%   NOTE: No checks are made regarding the correctness of the model. Use
%         checkModelStruct to identify problems in the model structure
%
%   Usage: exportToTabDelimited(model,path)
%
%   Cheewin Kittikunapong, 2019-04-02

if nargin<2
    path='./';
end

%If the folder doesn't exist then create it
if ~exist(path,'dir')
    mkdir(path);
end

%Remove the files if they already exist
if exist(fullfile(path,'excelRxns.txt'),'file')
    delete(fullfile(path,'excelRxns.txt'));
end
if exist(fullfile(path,'excelMets.txt'),'file')
    delete(fullfile(path,'excelMets.txt'));
end
if exist(fullfile(path,'excelGenes.txt'),'file')
    delete(fullfile(path,'excelGenes.txt'));
end
if exist(fullfile(path,'excelModel.txt'),'file')
    delete(fullfile(path,'excelModel.txt'));
end
if exist(fullfile(path,'excelComps.txt'),'file')
    delete(fullfile(path,'excelComps.txt'));
end

%Construct equations
model.equations=constructEquations(model,model.rxns,true);

%Open for printing the rxn sheet
rxnFile=fopen(fullfile(path,'excelRxns.txt'),'wt');

%Print header
fprintf(rxnFile,'#\tID\tNAME\tEQUATION\tEC-NUMBER\tGENE ASSOCIATION\tLOWER BOUND\tUPPER BOUND\tOBJECTIVE\tCOMPARTMENT\tMIRIAM\tSUBSYSTEM\tREPLACEMENT ID\tNOTE\tREFERENCE\tCONFIDENCE SCORE\n');

%Loop through the reactions
for i=1:numel(model.rxns)
    fprintf(rxnFile,['\t' model.rxns{i} '\t' model.rxnNames{i} '\t' model.equations{i} '\t']);
    
    if isfield(model,'eccodes')
        fprintf(rxnFile,[model.eccodes{i} '\t']);
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'grRules')
        fprintf(rxnFile,[model.grRules{i} '\t']);
    else
        fprintf(rxnFile,'\t');
    end
    
    %Print bounds and objectives
    fprintf(rxnFile,[num2str(model.lb(i)) '\t' num2str(model.ub(i)) '\t']);
    
    if model.c(i)~=0
        fprintf(rxnFile,[num2str(model.c(i)) '\t' ]);
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'rxnComps')
        fprintf(rxnFile,[model.comps{model.rxnComps(i)} '\t']);
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'rxnMiriams')
        if ~isempty(model.rxnMiriams{i})
            toPrint=[];
            for j=1:numel(model.rxnMiriams{i}.name)
                toPrint=[toPrint strtrim(model.rxnMiriams{i}.name{j}) '/' strtrim(model.rxnMiriams{i}.value{j}) ';'];
            end
            fprintf(rxnFile,[toPrint(1:end-1) '\t']);
        else
            fprintf(rxnFile,'\t');
        end
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'subSystems')
        if ~isempty(model.subSystems{i})
            fprintf(rxnFile,[strjoin(model.subSystems{i,1},';') '\t']);
        else
            fprintf(rxnFile,'\t');
        end
    end
    
    %Print replacement IDs
    fprintf(rxnFile,'\t');
    
    if isfield(model,'rxnNotes')
        fprintf(rxnFile,[model.rxnNotes{i} '\t']);
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'rxnReferences')
        fprintf(rxnFile,[model.rxnReferences{i} '\t']);
    else
        fprintf(rxnFile,'\t');
    end
    
    if isfield(model,'rxnConfidenceScores')
        fprintf(rxnFile,[num2str(model.rxnConfidenceScores(i)) '\t' ]);
    else
        fprintf(rxnFile,'\t');
    end
    
    fprintf(rxnFile,'\n');
end

fclose(rxnFile);

%Open for printing the metabolites sheet
metFile=fopen(fullfile(path,'excelMets.txt'),'wt');

%Print header
fprintf(metFile,'#\tID\tNAME\tUNCONSTRAINED\tMIRIAM\tCOMPOSITION\tInChI\tCOMPARTMENT\tREPLACEMENT ID\tMETS FIELD\tCHARGE\n');

%Loop through the metabolites
for i=1:numel(model.mets)
    fprintf(metFile,['\t' model.metNames{i} '[' model.comps{model.metComps(i)} ']\t' model.metNames{i} '\t']);
    
    if isfield(model,'unconstrained')
        if model.unconstrained(i)~=0
            fprintf(metFile,'TRUE\t');
        else
            fprintf(metFile,'\t');
        end
    else
        fprintf(metFile,'\t');
    end
    
    if isfield(model,'metMiriams')
        if ~isempty(model.metMiriams{i})
            toPrint=[];
            for j=1:numel(model.metMiriams{i}.name)
                toPrint=[toPrint strtrim(model.metMiriams{i}.name{j}) '/' strtrim(model.metMiriams{i}.value{j}) ';'];
            end
            fprintf(rxnFile,[toPrint(1:end-1) '\t']);
        else
            fprintf(metFile,'\t');
        end
    else
        fprintf(metFile,'\t');
    end
    
    if isfield(model,'metFormulas')
        fprintf(metFile,[model.metFormulas{i} '\t']);
    else
        fprintf(metFile,'\t');
    end
    
    if isfield(model,'inchis')
        fprintf(metFile,[model.inchis{i} '\t']);
    else
        fprintf(metFile,'\t');
    end
    
    fprintf(metFile,[model.comps{model.metComps(i)} '\t']);
    
    %There can be no replacement IDs in the structure, but it has to be
    %something to give working met IDs.
    fprintf(metFile,['m' int2str(i) '\t']);
    
    %Print the model.mets field. The reason for not putting this as
    %replacement ID is that it's not guaranteed to be a valid SBML id.
    fprintf(metFile,[model.mets{i} '\t']);
    
    if isfield(model,'metCharges')
        fprintf(metFile,[num2str(model.metCharges(i)) '\t']);
    else
        fprintf(metFile,'\t');
    end
    
    fprintf(metFile,'\n');
end

fclose(metFile);

if isfield(model,'genes')
    %Open for printing the genes sheet
    geneFile=fopen(fullfile(path,'excelGenes.txt'),'wt');
    
    %Print header
    fprintf(geneFile,'#\tNAME\tMIRIAM\tSHORT NAME\tCOMPARTMENT\n');
    
    %Loop through the genes
    for i=1:numel(model.genes)
        fprintf(geneFile,['\t' model.genes{i} '\t']);
        
        if isfield(model,'geneMiriams')
            if ~isempty(model.geneMiriams{i})
                toPrint=[];
                for j=1:numel(model.geneMiriams{i}.name)
                    toPrint=[toPrint strtrim(model.geneMiriams{i}.name{j}) '/' strtrim(model.geneMiriams{i}.value{j}) ';'];
                end
                fprintf(geneFile,[toPrint(1:end-1) '\t']);
            else
                fprintf(geneFile,'\t');
            end
        else
            fprintf(geneFile,'\t');
        end
        
        if isfield(model,'geneShortNames')
            fprintf(geneFile,[model.geneShortNames{i} '\t']);
        else
            fprintf(geneFile,'\t');
        end
        
        if isfield(model,'geneComps')
            fprintf(geneFile,[model.comps{model.geneComps(i)} '\t']);
        else
            fprintf(geneFile,'\t');
        end
        
        fprintf(geneFile,'\n');
    end
    fclose(geneFile);
end

if isfield(model,'id')
    %Open for printing the model sheet
    modelFile=fopen(fullfile(path,'excelModel.txt'),'wt');
    
    %Print header
    fprintf(geneFile,'#\tID\tDESCRIPTION\tDEFAULT LOWER\tDEFAULT UPPER\tCONTACT GIVEN NAME\tCONTACT FAMILY NAME\tCONTACT EMAIL\tORGANIZATION\tTAXONOMY\tNOTES\n');
    
    %Print model ID and name. It is assumed that the default lower/upper
    %bound correspond to min/max of the bounds
    toPrint=['\t' model.id '\t' model.description '\t'];
    if isfield(model,'annotation')
        if isfield(model.annotation,'defaultLB')
            toPrint=[toPrint num2str(model.annotation.defaultLB) '\t'];
        else
            toPrint=[toPrint num2str(min(model.lb)) '\t'];
        end
        if isfield(model.annotation,'defaultUB')
            toPrint=[toPrint num2str(model.annotation.defaultUB) '\t'];
        else
            toPrint=[toPrint num2str(max(model.ub)) '\t'];
        end
        if isfield(model.annotation,'givenName')
            toPrint=[toPrint model.annotation.givenName '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model.annotation,'familyName')
            toPrint=[toPrint model.annotation.familyName '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model.annotation,'email')
            toPrint=[toPrint model.annotation.email '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model.annotation,'organization')
            toPrint=[toPrint model.annotation.organization '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model.annotation,'taxonomy')
            toPrint=[toPrint model.annotation.taxonomy '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model.annotation,'note')
            toPrint=[toPrint model.annotation.note '\t'];
        else
            toPrint=[toPrint '\t'];
        end
    else
        toPrint=[toPrint num2str(min(model.lb)) '\t' num2str(max(model.ub)) '\tRasmus\tAgren\trasmus.agren@scilifelab.se\tChalmers University of Technology\t\t\n'];
    end
    fprintf(modelFile,toPrint);
    fclose(modelFile);
end

if isfield(model,'comps')
    %Open for printing the model sheet
    compsFile=fopen(fullfile(path,'excelComps.txt'),'wt');
    
    %Print header
    fprintf(compsFile,'#\tABBREVIATION\tNAME\tINSIDE\tMIRIAM\n');
    
    for i=1:numel(model.comps)
        toPrint=['\t' model.comps{i} '\t' model.compNames{i} '\t'];
        if isfield(model,'compOutside')
            toPrint=[toPrint model.compOutside{i} '\t'];
        else
            toPrint=[toPrint '\t'];
        end
        if isfield(model,'compMiriams')
            if ~isempty(model.compMiriams{i})
                for j=1:numel(model.compMiriams{i}.name)
                    toPrint=[toPrint strtrim(model.compMiriams{i}.name{j}) '/' strtrim(model.compMiriams{i}.value{j}) ';'];
                end
                toPrint(end)=[];
                toPrint=[toPrint '\t'];
            else
                toPrint=[toPrint '\t'];
            end
        else
            toPrint=[toPrint '\t'];
        end
        toPrint=[toPrint '\n'];
        fprintf(compsFile,toPrint);
    end
    fclose(compsFile);
end
end
