function model=importExcelModel(fileName,removeExcMets,printWarnings,ignoreErrors)
% importExcelModel
%   Imports a constraint-based model from a Excel file
%
%   fileName      a Microsoft Excel file to import
%   removeExcMets true if exchange metabolites should be removed. This is
%                 needed to be able to run simulations, but it could also 
%                 be done using simplifyModel at a later stage (opt,
%                 default true)
%   printWarnings true if warnings should be printed (opt, default true)
%   ignoreErrors  true if errors should be ignored. See below for details
%                 (opt, default false)
%
%   model
%       description      description of model contents
%       id               model ID
%       rxns             reaction ids
%       mets             metabolite ids
%       S                stoichiometric matrix
%       lb               lower bounds
%       ub               upper bounds
%       rev              reversibility vector
%       c                objective coefficients
%       b                equality constraints for the metabolite equations
%       comps            compartment ids
%       compNames        compartment names
%       compOutside      the id (as in comps) for the compartment
%                        surrounding each of the compartments
%       compMiriams      structure with MIRIAM information about the
%                        compartments
%       rxnNames         reaction description
%       rxnComps         compartments for reactions
%       grRules          reaction to gene rules in text form
%       rxnGeneMat       reaction-to-gene mapping in sparse matrix form
%       subSystems       subsystem name for each reaction
%       eccodes          EC-codes for the reactions
%       rxnMiriams       structure with MIRIAM information about the reactions
%       rxnNotes         reaction notes
%       rxnReferences	 reaction references
%       confidenceScores reaction confidence scores
%       genes            list of all genes
%       geneComps        compartments for reactions
%       geneMiriams      structure with MIRIAM information about the genes
%       metNames         metabolite description
%       metComps         compartments for metabolites
%       inchis           InChI-codes for metabolites
%       metFormulas      metabolite chemical formula
%       metMiriams       structure with MIRIAM information about the metabolites
%       metCharge        metabolite charge
%       unconstrained    true if the metabolite is an exchange metabolite
%
%   Loads models in the RAVEN Toolbox Excel format. A number of consistency
%   checks are performed in order to ensure that the model is valid. These
%   can be ignored by putting ignoreErrors to true. However, this is highly
%   advised against, as it can result in errors in simulations or other
%   functionalities. The RAVEN Toolbox is made to function only on consistent
%   models, and the only checks performed are when the model is imported.
%
%   NOTE: Most errors are checked for by checkModelStruct, but some
%   are checked for in this function as well. Those are ones which relate
%   to missing model elements and so on, and which would make it impossible
%   to construct the model structure. Those errors cannot be ignored by
%   setting ignoreErrors to true.
%
%   Usage: model=importExcelModel(fileName,removeExcMets,printWarnings,ignoreErrors)
%
%   Rasmus Agren, 2014-01-06
%   Simonas Marcisauskas, 2016-11-01 - added support for rxnNotes,
%   rxnReferences, confidenceScores and metCharge
%

%Adds the required classes to the Java path
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));
% Adding escape characters, if some parent folders contain spaces or
% exclamation marks (for Unix systems). For Windows, all the parent folders
% are just put between the double quotation brackets
if isunix
    ravenPath = regexprep(ravenPath,'\ ','\\ ');
    ravenPath = regexprep(ravenPath,'\!','\\!');
elseif ispc
    for i=1:(length(strfind(ravenPath,'\')))
        if i==1
            ravenPath = regexprep(ravenPath,'\\','\\"',i);
        elseif i==length(strfind(ravenPath,'\'))
            ravenPath = regexprep(ravenPath,'\\','"\\',i);    
        else
            ravenPath = regexprep(ravenPath,'\\','"\\"',i);
        end
    end
end

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
import java.io.FileInputStream;

if nargin<2
    removeExcMets=true;
end
if nargin<3
    printWarnings=true;
end
if nargin<4
    ignoreErrors=false;
end

%This is to match the order of the fields to those you get from importing
%from SBML
model=[];
model.id=[];
model.description=[];
model.annotation=[];
%Default bounds if not defined
model.annotation.defaultLB=-1000;
model.annotation.defaultUB=1000;
model.rxns={};
model.mets={};
model.S=[];
model.lb=[];
model.ub=[];
model.rev=[];
model.c=[];
model.b=[];
model.comps={};
model.compNames={};
model.compOutside={};
model.compMiriams={};
model.rxnNames={};
model.rxnComps={}; %Will be double later
model.grRules={};
model.rxnGeneMat=[];
model.subSystems={};
model.eccodes={};
model.rxnMiriams={};
model.rxnNotes={};
model.rxnReferences={};
model.confidenceScores={};
model.genes={};
model.geneComps={}; %Will be double later
model.geneMiriams={};
model.metNames={};
model.metComps=[];
model.inchis={};
model.metFormulas={};
model.metMiriams={};
model.metCharge={}; %Will be double later
model.unconstrained=[];
    
%Check if the file exists
if ~exist(fileName,'file')
    dispEM('The Excel file could not be found');
end

%Opens the workbook. The input stream is needed since it will otherwise
%keep the file open
is=FileInputStream(getFullPath(fileName));
workbook=WorkbookFactory.create(is);
is.close();

[raw, flag]=loadSheet(workbook,'MODEL');

if flag<0
    if printWarnings==true
        dispEM('Could not load the MODEL sheet',false);
    end
    model.id='UNKNOWN';
    model.description='No model details available';
else
    raw=cleanImported(raw);

    %It is assumed that the first line is labels and that the second one is
    %info
    allLabels={'ID';'DESCRIPTION';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'TAXONOMY';'NOTES'};

    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'MODELID','ID');
    raw(1,:)=strrep(raw(1,:),'MODELNAME','DESCRIPTION');

    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
                if any(raw{I(i),2})
                    model.id=toStr(raw{2,I(i)}); %Should be string already
                else
                    dispEM('No model ID supplied');
                end
            case 2
                if any(raw{2,I(i)})
                    model.description=toStr(raw{2,I(i)}); %Should be string already
                else
                    dispEM('No model name supplied');
                end
            case 3
                if ~isempty(raw{2,I(i)})
                    try
                        model.annotation.defaultLB=toDouble(raw{2,I(i)},NaN);
                    catch
                        dispEM('DEFAULT LOWER must be numeric');
                    end
                else
                    if printWarnings==true
                        fprintf('NOTE: DEFAULT LOWER not supplied. Uses -1000\n');
                    end
                    model.annotation.defaultLB=-1000;
                end
            case 4
                if ~isempty(raw{2,I(i)})
                    try
                        model.annotation.defaultUB=toDouble(raw{2,I(i)},NaN);
                    catch
                        dispEM('DEFAULT UPPER must be numeric');
                    end
                else
                    if printWarnings==true
                        fprintf('NOTE: DEFAULT UPPER not supplied. Uses 1000\n');
                    end
                    model.annotation.defaultUB=1000;
                end
            case 5
                if any(raw{2,I(i)})
                    model.annotation.givenName=toStr(raw{2,I(i)}); %Should be string already 
                end
            case 6
                if any(raw{2,I(i)})
                    model.annotation.familyName=toStr(raw{2,I(i)}); %Should be string already 
                end
            case 7
                if any(raw{2,I(i)})
                    model.annotation.email=toStr(raw{2,I(i)}); %Should be string already 
                end
            case 8
                if any(raw{2,I(i)})
                    model.annotation.organization=toStr(raw{2,I(i)}); %Should be string already 
                end    
            case 9
                if any(raw{2,I(i)})
                    model.annotation.taxonomy=toStr(raw{2,I(i)}); %Should be string already 
                end 
            case 10
                if any(raw{2,I(i)})
                    model.annotation.note=toStr(raw{2,I(i)}); %Should be string already 
                end     
        end  
    end
end

%Get compartment information
[raw, flag]=loadSheet(workbook,'COMPS');

if flag<0
    if printWarnings==true
        dispEM('Could not load the COMPS sheet. All elements will be assigned to a compartment "s" for "System"',false);
    end
    model.comps={'s'};
    model.compNames={'System'};
else
    raw=cleanImported(raw);

    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'COMPABBREV','ABBREVIATION');
    raw(1,:)=strrep(raw(1,:),'COMPNAME','NAME');

    allLabels={'ABBREVIATION';'NAME';'INSIDE';'MIRIAM'};

    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
                model.comps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 2
                model.compNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                model.compOutside=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 4
                model.compMiriams=raw(2:end,I(i));
        end
    end

    %Check that necessary fields are loaded
    if isempty(model.comps)
        dispEM('There must be a column named ABBREVIATION in the COMPS sheet');   
    end
    if isempty(model.compNames)
        model.compNames=model.comps; 
        if printWarnings==true
            dispEM('There is no column named NAME in the COMPS sheet. ABBREVIATION will be used as name',false);
        end
    end

    model.compMiriams=parseMiriam(model.compMiriams);
end

%Get all the genes and info about them
[raw, flag]=loadSheet(workbook,'GENES');

if flag<0
    if printWarnings==true
        dispEM('There is no spreadsheet named GENES',false)
    end
else
    raw=cleanImported(raw);

    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'GENE NAME','NAME');

    allLabels={'NAME';'MIRIAM';'COMPARTMENT'};
    
    %Loop through the labels
    [I, J]=ismember(upper(raw(1,:)),allLabels);
    I=find(I);
    foundGenes=false;
    for i=1:numel(I)
        switch J(I(i))
            case 1
                model.genes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
                foundGenes=true;
            case 2
                model.geneMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                model.geneComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        end
    end

    if foundGenes==false
        dispEM('There must be a column named NAME in the GENES sheet');
    end
        
    %Its ok if all of them are empty
    if all(cellfun(@isempty,model.geneComps))
        model.geneComps=[];
    end

    %Check that geneName contain only strings and no empty strings
    if ~iscellstr(model.genes)
        dispEM('All gene names have to be strings');
    else
        if any(strcmp('',model.genes))
            dispEM('There can be no empty strings in gene names');
        end
    end

    %Check that geneComp contains only strings and no empty string    
    if ~isempty(model.geneComps)
        if ~iscellstr(model.geneComps)
            dispEM('All gene compartments have to be strings');
        else
            if any(strcmp('',model.geneComps))
                dispEM('There can be no empty strings in gene compartments');
            end
        end
        [I, model.geneComps]=ismember(model.geneComps,model.comps);
        dispEM('The following genes have compartment abbreviations which could not be found:',true,model.genes(~I));
    end
end

model.geneMiriams=parseMiriam(model.geneMiriams);

%Loads the reaction data
[raw, flag]=loadSheet(workbook,'RXNS');

if flag<0
    dispEM('Could not load the RXNS sheet');
end

raw=cleanImported(raw);

%Map to new captions
raw(1,:)=upper(raw(1,:));
raw(1,:)=strrep(raw(1,:),'RXNID','ID');

allLabels={'ID';'NAME';'EQUATION';'EC-NUMBER';'GENE ASSOCIATION';'LOWER BOUND';'UPPER BOUND';'OBJECTIVE';'COMPARTMENT';'SUBSYSTEM';'REPLACEMENT ID';'MIRIAM';'NOTE';'REFERENCE';'CONFIDENCE SCORE'};
equations={};
reactionReplacement={};

%Loop through the labels
[I, J]=ismember(raw(1,:),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
           model.rxns=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 2
           model.rxnNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 3
           equations=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 4
           model.eccodes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 5
           model.grRules=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 6
           try
               model.lb=cellfun(@(x) toDouble(x,NaN),raw(2:end,I(i)));
           catch
               dispEM('The lower bounds must be numerical values'); 
           end
        case 7
           try
               model.ub=cellfun(@(x) toDouble(x,NaN),raw(2:end,I(i)));
           catch
               dispEM('The upper bounds must be numerical values'); 
           end
        case 8
           try
               model.c=cellfun(@(x) toDouble(x,0),raw(2:end,I(i)));
           catch
               dispEM('The objective coefficients must be numerical values'); 
           end
        case 9
        	model.rxnComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 10
        	model.subSystems=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 11
        	reactionReplacement=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);    
        case 12
            model.rxnMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 13
            model.rxnNotes=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 14
            model.rxnReferences=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        case 15
            model.confidenceScores=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
    end
end

%Check that all necessary reaction info has been loaded
if isempty(equations)
     dispEM('There must be a column named EQUATION in the RXNS sheet');   
end
if isempty(model.rxns)
     if printWarnings==true
        dispEM('There is no column named ID in the RXNS sheet. The reactions will be named as "r1", "r2"...',false);
     end
     I=num2cell((1:numel(equations))');
     model.rxns=strcat('r',cellfun(@num2str,I,'UniformOutput',false));
end

%Check if some other stuff is loaded and populate with default values
%otherwise
if isempty(model.rxnNames)
   model.rxnNames=cell(numel(model.rxns),1);
   model.rxnNames(:)={''};
   if printWarnings==true
       dispEM('There is no column named NAME in the RXNS sheet. Empty strings will be used as reaction names',false);
   end
end
if isempty(model.lb)
   %This is not set here since the reversibility isn't known yet 
   model.lb=nan(numel(model.rxns),1);
   if printWarnings==true
        dispEM('There is no column named LOWER BOUND in the RXNS sheet. Default bounds will be used',false);
   end
end
if isempty(model.ub)
   model.ub=nan(numel(model.rxns),1);
   if printWarnings==true
        dispEM('There is no column named UPPER BOUND in the RXNS sheet. Default bounds will be used',false);
   end
end
if isempty(model.c)
   model.c=zeros(numel(model.rxns),1);
   if printWarnings==true
        dispEM('There is no column named OBJECTIVE in the RXNS sheet',false);
   end
end

%Either all reactions must have a compartment string or none of them.
%Check if it's only empty and if so return it to []
if ~isempty(model.rxnComps)
    if all(cellfun(@isempty,model.rxnComps))
        model.rxnComps=[];
    end
end

%Construct the rxnMiriams structure
model.rxnMiriams=parseMiriam(model.rxnMiriams);

%Replace the reaction IDs for those IDs that have a corresponding 
%replacement name.
I=cellfun(@any,reactionReplacement);
model.rxns(I)=reactionReplacement(I);

%Check that there are no empty strings in reactionIDs or equations
if any(strcmp('',model.rxns))
    dispEM('There are empty reaction IDs'); 
end

if any(strcmp('',equations))
    dispEM('There are empty equations'); 
end

if ~isempty(model.rxnComps)
    if any(strcmp('',model.rxnComps))
        dispEM('Either all reactions must have an associated compartment string or none of them'); 
    end
end
    
%Check gene association for each reaction and populate rxnGeneMat
if ~isempty(model.genes)
    model.rxnGeneMat=zeros(numel(model.rxns),numel(model.genes));
end
if ~isempty(model.grRules)
    for i=1:length(model.rxns)    
       %Check that all gene associations have a match in the gene list
       if ~isempty(model.grRules{i})
           indexes=strfind(model.grRules{i},':'); %Genes are separated by ":" for AND and ";" for OR
           indexes=unique([indexes strfind(model.grRules{i},';')]);
           if isempty(indexes)
               %See if you have a match
               I=find(strcmp(model.grRules{i},model.genes));
               if isempty(I)
                   dispEM(['The gene association in reaction ' model.rxns{i} ' (' model.grRules{i} ') is not present in the gene list']);
               end
               model.rxnGeneMat(i,I)=1;
           else
               temp=[0 indexes numel(model.grRules{i})+1];
               for j=1:numel(indexes)+1;
                   %The reaction has several associated genes
                   geneName=model.grRules{i}(temp(j)+1:temp(j+1)-1);
                   I=find(strcmp(geneName,model.genes));
                   if isempty(I)
                       dispEM(['The gene association in reaction ' model.rxns{i} ' (' geneName ') is not present in the gene list']);
                   end
                   model.rxnGeneMat(i,I)=1;
               end
           end
            %In order to adhere to the COBRA standards it should be like
            %this:
            %-If only one gene then no parentheses
            %-If only "and" or only "or" there should only be one set of
            %parentheses
            %-If both "and" and "or", then split on "or". This is not
            %complete, but it's the type of relationship supported by the
            %Excel formulation
            aSign=strfind(model.grRules{i},':');
            oSign=strfind(model.grRules{i},';');
            if isempty(aSign) && isempty(oSign)
                model.grRules{i}=model.grRules{i};
            else
                if isempty(aSign)
                    model.grRules{i}=['(' strrep(model.grRules{i},';',' or ') ')'];
                else
                    if isempty(oSign)
                        model.grRules{i}=['(' strrep(model.grRules{i},':',' and ') ')'];
                    else
                        model.grRules{i}=['((' strrep(model.grRules{i},';',') or (') '))'];
                        model.grRules{i}=strrep(model.grRules{i},':',' and ');
                    end
                end
            end
       end
    end
end
model.rxnGeneMat=sparse(model.rxnGeneMat);

%Check that the compartment for each reaction can be found
if ~isempty(model.rxnComps)
    [I, model.rxnComps]=ismember(model.rxnComps,model.comps);
    dispEM('The following reactions have compartment abbreviations which could not be found:',true,model.rxns(~I));
end

%Get all the metabolites and info about them
[raw, flag]=loadSheet(workbook,'METS');

if flag<0
    if printWarnings==true
        dispEM('There is no spreadsheet named METS. The metabolites will be named "m1", "m2"... and assigned to the first compartment',false);
    end
    %Parse the equations to find out how many metabolites there are
    metsForParsing=parseRxnEqu(equations);
    I=num2cell((1:numel(metsForParsing))');
    model.mets=strcat('m',cellfun(@num2str,I,'UniformOutput',false));
    model.metComps=ones(numel(model.mets),1);
    model.unconstrained=zeros(numel(model.mets),1);
    model.metNames=metsForParsing;
else
    raw=cleanImported(raw);			

    %Map to new captions
    raw(1,:)=upper(raw(1,:));
    raw(1,:)=strrep(raw(1,:),'METID','ID');
    raw(1,:)=strrep(raw(1,:),'METNAME','NAME');

    allLabels={'ID';'NAME';'UNCONSTRAINED';'MIRIAM';'COMPOSITION';'INCHI';'COMPARTMENT';'REPLACEMENT ID';'CHARGE'};

    %Load the metabolite information
    metReplacement={};

    %Loop through the labels
    [I, J]=ismember(raw(1,:),allLabels);
    I=find(I);
    for i=1:numel(I)
        switch J(I(i))
            case 1
            	model.mets=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 2
            	model.metNames=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 3
            	model.unconstrained=cellfun(@boolToDouble,raw(2:end,I(i)));
               
            	%NaN is returned if the values couldn't be parsed
            	dispEM('The UNCONSTRAINED property for the following metabolites must be "true"/"false", 1/0, TRUE/FALSE or not set:',true,model.mets(isnan(model.unconstrained)));
            case 4
            	model.metMiriams=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 5
            	model.metFormulas=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 6
            	model.inchis=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 7
                model.metComps=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
                %Check that all metabolites have compartments defined
                if any(strcmp('',model.metComps))
                    dispEM('All metabolites must have an associated compartment string'); 
                end
            case 8
            	metReplacement=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
            case 9
            	model.metCharge=cellfun(@toStr,raw(2:end,I(i)),'UniformOutput',false);
        end
    end

    %Check that necessary fields are loaded (METID)
    if isempty(model.mets)
         dispEM('There must be a column named ID in the METS sheet');   
    end

    %Check that some other stuff is loaded and use default values otherwise
    if isempty(model.metNames)
       model.metNames=cell(numel(model.mets),1);
       if printWarnings==true
            dispEM('There is no column named NAME in the METS sheet. ID will be used as name',false);
       end
    end
    if isempty(model.unconstrained)
       model.unconstrained=zeros(numel(model.mets),1);
       if printWarnings==true
            dispEM('There is no column named UNCONSTRAINED in the METS sheet. All metabolites will be constrained',false);
       end
    end

    if isempty(model.metComps)
       model.metComps=cell(numel(model.mets),1);
       model.metComps(:)=model.comps(1);
       if printWarnings==true
            dispEM('There is no column named COMPARTMENT in the METS sheet. All metabolites will be assigned to the first compartment in COMPS. Note that RAVEN makes extensive use of metabolite names and compartments. Some features will therefore not function correctly if metabolite compartments are not correctly assigned',false);
       end
    end

    %The composition should be loaded from InChIs when available
    I=find(~cellfun(@isempty,model.inchis));
    for i=1:numel(I)
        S=regexp(model.inchis(I(i)),'/','split');
        S=S{1};
        if numel(S)>=2
            %Don't copy if it doesn't look good
            model.metFormulas(I(i))=S(2);
        end
    end

    %Check that the compartment for each metabolite can be found. Also convert
    %from id to index
    [I, model.metComps]=ismember(model.metComps,model.comps);
    dispEM('The following metabolites have compartment abbreviations which could not be found:',true,model.mets(~I));
    
    %Check that the model.mets vector is unique. The problem is that
    %the checkModelStruct cannot check for that since only the metReplacements
    %(if used) end up in the model structure
    I=false(numel(model.mets),1);
    [J, K]=unique(model.mets);
    if numel(J)~=numel(model.mets)
        L=1:numel(model.mets);
        L(K)=[];
        I(L)=true;
    end
    dispEM('The following metabolites are duplicates:',~ignoreErrors,model.mets(I));

    %Check that there are no metabolite IDs which are numbers. This would
    %give errors when parsing the equations
    I=cellfun(@str2double,model.mets);
    dispEM('The following metabolites have names which cannot be distinguished from numbers:',~ignoreErrors,model.mets(~isnan(I)));
    I=cellfun(@str2double,metReplacement);
    dispEM('The following metabolites have names which cannot be distinguished from numbers:',~ignoreErrors,metReplacement(~isnan(I)));
    
    %Replace the metabolite IDs for those IDs that have a corresponding 
    %replacement metabolite. This is not used for matching, but will be checked
    %for consistency with SBML naming conventions
    metsForParsing=model.mets; %This is because the equations are written with this 
    I=cellfun(@any,metReplacement);
    model.mets(I)=metReplacement(I);
    
    %If the metabolite name isn't set, replace it with the metabolite id
    I=~cellfun(@any,model.metNames);
    model.metNames(I)=model.mets(I);

    %Construct the metMiriams structure
    model.metMiriams=parseMiriam(model.metMiriams);
    
    %Either all metabolites have charge or none of them.
    %Check if it's only empty and if so return it to []
    if ~isempty(model.metCharge)
        if all(cellfun(@isempty,model.metCharge))
            model.metCharge=[];
        end
    end
    if ~isempty(model.metCharge)
        if any(strcmp('',model.metCharge))
            dispEM('Either all metabolites have charge information or none of them'); 
        end
    end
    if ~isempty(model.metCharge)
        model.metCharge=str2double(model.metCharge);
    end
end

%Everything seems fine with the metabolite IDs, compartments, genes, and
%reactions

%Parse the equations
[model.S, mets, badRxns, model.rev]=constructS(equations,metsForParsing,model.rxns);
model.rev=model.rev*1; %Typecast to double

%Add default constraints
model.lb(isnan(model.lb))=model.annotation.defaultLB.*model.rev(isnan(model.lb));
model.ub(isnan(model.ub))=model.annotation.defaultUB;

%Reorder the S matrix so that it fits with the metabolite list in the
%structure
[crap, I]=ismember(mets,metsForParsing);
model.S=model.S(I,:);

%Print warnings about the reactions which contain the same metabolite as
%both reactants and products
dispEM('The following reactions have metabolites which are present more than once. Only the net reactions will be exported:',false,model.rxns(badRxns));

model.b=zeros(numel(model.mets),1);

%Remove unused fields
if isempty(model.compOutside)
    model=rmfield(model,'compOutside');
end
if isempty(model.compMiriams)
    model=rmfield(model,'compMiriams');
end
if isempty(model.rxnComps)
    model=rmfield(model,'rxnComps');
end
if isempty(model.grRules)
    model=rmfield(model,'grRules');
end
if isempty(model.rxnGeneMat)
    model=rmfield(model,'rxnGeneMat');
end
if isempty(model.subSystems)
    model=rmfield(model,'subSystems');
end
if isempty(model.eccodes)
    model=rmfield(model,'eccodes');
end
if isempty(model.rxnMiriams)
	model=rmfield(model,'rxnMiriams');
end
if isempty(model.rxnNotes)
	model=rmfield(model,'rxnNotes');
end
if isempty(model.rxnReferences)
	model=rmfield(model,'rxnReferences');
end
if isempty(model.confidenceScores)
	model=rmfield(model,'confidenceScores');
end
if isempty(model.genes)
    model=rmfield(model,'genes');
end
if isempty(model.geneComps)
    model=rmfield(model,'geneComps');
end
if isempty(model.geneMiriams)
    model=rmfield(model,'geneMiriams');
end
if isempty(model.inchis)
    model=rmfield(model,'inchis');
end
if isempty(model.metFormulas)
    model=rmfield(model,'metFormulas');
end
if isempty(model.metMiriams)
    model=rmfield(model,'metMiriams');
end
if isempty(model.metCharge)
    model=rmfield(model,'metCharge');
end

%The model structure has now been reconstructed but it can still contain
%many types of errors. The checkModelConsistency function is used to make
%sure that naming and mapping of stuff looks good
checkModelStruct(model,~ignoreErrors);

if removeExcMets==true
    model=simplifyModel(model);
end
end

%Cleans up the structure that is imported from using xlsread
function raw=cleanImported(raw)
    %First check that it's a cell array. If a sheet is completely empty,
    %then raw=NaN
    if iscell(raw)
        %Clear cells which contain only white spaces or NaN. This could
        %happen if you accidentally inserted a space for example. I don't
        %know how NaN could occur after switching to Apache POI, but I
        %clear it to be sure
        whites=cellfun(@wrapperWS,raw);
        raw(whites)={[]};
        
        %Remove columns that don't have string headers. If you cut and paste
        %a lot in the sheet there tends to be columns that are empty
        I=cellfun(@isstr,raw(1,:));
        raw=raw(:,I);
        
        %Find the rows that are not commented. This corresponds to the
        %first row and the ones which are empty in the first column
        
        keepers=cellfun(@isempty,raw(:,1));
        keepers(1)=true;
        raw=raw(keepers,:);

        %Check if there are any rows that are all empty. This could happen if
        %it reads too far or if the user has inserted them for layout reasons.
        %Remove any such rows.
        I=~all(cellfun(@isempty,raw),2);
        raw=raw(I,:);
    else
        raw={[]};
    end
    
    %Checks if something is all white spaces or NaN
    function I=wrapperWS(A)
        if isnan(A)
            I=true;
        else
            %isstrprop gives an error if boolean
            if islogical(A)
                I=false;
            else
                %I have to make this check, since otherwise it will pick up
                %on doubles for which the ASCII representation is a white
                %space character
                if isnumeric(A)
                    I=false;
                else
                    I=all(isstrprop(A,'wspace'));
                end
            end
        end
    end
end
function miriamStruct=parseMiriam(strings,miriamStruct)
%Gets the names and values of Miriam-string. Nothing fancy at all, just to
%prevent using the same code for metabolites, genes, and reactions. The
%function also allows for supplying a miriamStruct and the info will then
%be added

if nargin<2
    miriamStruct=cell(numel(strings),1);
end
for i=1:numel(strings)
    if any(strings{i})
        %A Miriam string can be several ids separated by ";". Each id is
        %"name(..:..):value"
        I=regexp(strings{i},';','split');
        if isfield(miriamStruct{i},'name')
            startIndex=numel(miriamStruct{i}.name);
            miriamStruct{i}.name=[miriamStruct{i}.name;cell(numel(I),1)];
            miriamStruct{i}.value=[miriamStruct{i}.value;cell(numel(I),1)];
        else
            startIndex=0;
            miriamStruct{i}.name=cell(numel(I),1);
            miriamStruct{i}.value=cell(numel(I),1);
        end
        
        for j=1:numel(I)
            index=max(strfind(I{j},':'));
            if any(index)
                miriamStruct{i}.name{startIndex+j}=I{j}(1:index-1);
                miriamStruct{i}.value{startIndex+j}=I{j}(index+1:end);
            else
                dispEM(['"' I{j} '" is not a valid MIRIAM string. The format must be "identifier:value"']);
            end
        end
    end
end
end

%Loads a sheet into a cell matrix using the Java library Apache POI.

%flag       -1 if the sheet couldn't be found
%           0 if everything worked well
function [raw, flag]=loadSheet(workbook, sheet)      
    flag=0;
    raw={};
    
    sh=workbook.getSheet(sheet);
    if isempty(sh)
        flag=-1;
        return;
    end
    
    lastRow=sh.getLastRowNum();
    wasEmpty=false(lastRow+1,1);
    raw=cell(lastRow+1,0); %Allocate space for the cell array. The number of columns isn't know yet, as it's saved row-by-row
    for i=0:lastRow
        row=sh.getRow(i);
        %Sometimes the last rows only contain formatting (or some other
        %weird Excel thing). Ignore such empty rows. Note the +1 to deal with that Matlab indexing
        %starts at 1
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
               %Then decide how to save it depending on the type
               switch c.getCellType()
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

%For converting a value to string. This is used instead of num2str because
%I want to convert empty cells to {''}.
function y=toStr(x)
    %x can be empty, numerical, string or boolean. It cannot be NaN.
    %Boolean values will be converted to '1'/'0'
    if isempty(x)
        y='';
    else
        y=num2str(x);
    end
end

%For converting to numeric. This is used instead of str2num because I want
%to be able to choose what empty values should be mapped to.
%
% default the value to use for empty input
function y=toDouble(x,default)
    if isempty(x) %Note that this catches '' as well
        y=default;
    else
        if isnumeric(x)
            y=x;
        else
            y=str2double(x);
            
            %This happens if the input couldn't be converted.
            %Note that the input itself cannot be NaN since it was fixed in
            %clean imported
            if isnan(y)
                dispEM(['Cannot convert the string "' x '" to double']);
            end
        end
    end
end

%For converting boolean (the UNCONSTRAINED field) to double (the
%model.unconstrained field)
function y=boolToDouble(x)
    if isempty(x)
        y=0;
        return;
    end
    if islogical(x)
        y=x*1; %Typecast to double
        return;
    end
    if isnumeric(x)
       if x~=0
           y=1;
           return;
       else
           y=0;
           return;
       end
    end
    if ischar(x)
       if strcmpi(x,'TRUE')
           y=1;
           return;
       end
       if strcmpi(x,'FALSE')
           y=0;
           return;
       end
    end
    y=NaN; %This means that the input couldn't be parsed
end