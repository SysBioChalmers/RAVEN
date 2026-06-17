function exportToExcelFormat(model,varargin)
% exportToExcelFormat  Export a model to the Microsoft Excel model format.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% fileName : char
%     file name of the Excel file. Only xlsx format is supported. This can
%     also be only a path, in which case the model is exported to a set of
%     tab-delimited text files via exportToTabDelimited. A dialog window
%     will open if fileName is empty.
% sortIds : logical
%     whether metabolites, reactions and genes should be sorted
%     alphabetically by their identifiers (default false).
%
% Examples
% --------
%     exportToExcelFormat(model, fileName, sortIds);

p=parseRAVENargs(varargin, {'fileName',[]; 'sortIds',false});
fileName=p.fileName; sortIds=p.sortIds;
if isempty(fileName)
    [fileName, pathName] = uiputfile('*.xlsx', 'Select file for model export',[model.id '.xlsx']);
    if fileName == 0
        error('You should provide a file location')
    else
        fileName = fullfile(pathName,fileName);
    end
end
fileName=char(fileName);
if sortIds==true
    model=sortIdentifiers(model);
end

checkModelStruct(model);

[~, A, B]=fileparts(fileName);

%If a path was used call on exportToTabDelimited instead
if ~any(A) || ~any(B)
    exportToTabDelimited(model,fileName);
    return;
end

if ~strcmpi(B,'.xlsx')
    EM='As of RAVEN version 1.9, only export to xlsx format is supported';
    dispEM(EM);
end

%Accumulate the sheets, then write them all to the .xlsx at the end.
sheets=struct('name',{},'header',{},'data',{},'widths',{},'isIntegers',{},'colFormats',{});

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

subsystems='';
if isfield(model,'subSystems')
    for i=1:numel(model.subSystems)
        if ~iscell(model.subSystems{i})
            model.subSystems{i} = {model.subSystems{i}};
        end    
        if ~isempty(model.subSystems{i,1})
            subsystems{i,1}=strjoin(model.subSystems{i,1},';');
        else
            subsystems{i,1}='';
        end
    end
    rxnSheet=[rxnSheet subsystems];
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
    rxnSheet=[rxnSheet num2cell(model.rxnConfidenceScores)];
else
    rxnSheet=[rxnSheet emptyColumn];
end

sheets=appendSheet(sheets,'RXNS',headers,rxnSheet,true,{});

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
    
    % Making sure that only these metFormulas are exported, which don't
    % have InChI strings
    if isfield(model,'metFormulas')
        if isfield(model,'inchis')
            if isempty(model.inchis{i})
                metSheet(i,6)=model.metFormulas(i);
            end
        else
            metSheet(i,6)=model.metFormulas(i);
        end
    end
    
    if isfield(model,'inchis')
        metSheet(i,7)=model.inchis(i);
    end
    
    if isfield(model,'metComps')
        metSheet(i,8)=model.comps(model.metComps(i));
    end
    
    metSheet(i,9)=model.mets(i);
    
    if isfield(model,'metCharges')
        metSheet{i,10}=model.metCharges(i);
    end
end

sheets=appendSheet(sheets,'METS',headers,metSheet,true,{});

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

sheets=appendSheet(sheets,'COMPS',headers,compSheet,true,{});

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
    
    sheets=appendSheet(sheets,'GENES',headers,geneSheet,true,{});
end

%Add the MODEL sheet

%Create the header row
headers={'#';'ID';'NAME';'TAXONOMY';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'NOTES'};

modelSheet=cell(1,numel(headers));

if ~isfield(model,'annotation')
    model.annotation = [];
end

if isfield(model,'id')
    modelSheet{1,2}=model.id;
else
    modelSheet{1,2}='blankID';
end
if isfield(model,'name')
    modelSheet{1,3}=model.name;
else
    modelSheet{1,3}='blankName';
end
if isfield(model.annotation,'taxonomy')
    modelSheet{1,4}=model.annotation.taxonomy;
end
if isfield(model.annotation,'defaultLB')
    modelSheet{1,5}=model.annotation.defaultLB;
end
if isfield(model.annotation,'defaultUB')
    modelSheet{1,6}=model.annotation.defaultUB;
end
if isfield(model.annotation,'givenName')
    modelSheet{1,7}=model.annotation.givenName;
end
if isfield(model.annotation,'familyName')
    modelSheet{1,8}=model.annotation.familyName;
end
if isfield(model.annotation,'email')
    modelSheet{1,9}=model.annotation.email;
end
if isfield(model.annotation,'organization')
    modelSheet{1,10}=model.annotation.organization;
end
if isfield(model.annotation,'note')
    modelSheet{1,11}=model.annotation.note;
end

sheets=appendSheet(sheets,'MODEL',headers,modelSheet,true,{});

%Add the ENZYMES and ENZRXNS sheets, containing the contents of the
%model.ec structure of enzyme-constrained (GECKO) models. The
%enzyme-reaction coupling (model.ec.rxnEnzMat) is written as the
%'enzyme:count' ENZYMES column of the ENZRXNS sheet.
if isfield(model,'ec') && isfield(model.ec,'enzymes')
    %ENZYMES sheet: one row per enzyme
    headers={'#';'ID';'GENE';'MW';'SEQUENCE';'CONC'};
    enzSheet=cell(numel(model.ec.enzymes),numel(headers));
    enzSheet(:,2)=model.ec.enzymes(:);
    enzSheet(:,3)=model.ec.genes(:);
    enzSheet(:,4)=numToCell(model.ec.mw);
    enzSheet(:,5)=model.ec.sequence(:);
    enzSheet(:,6)=numToCell(model.ec.concs);
    %MW (column 4) shown without decimals, CONC (column 6) with 5 decimals
    enzColFormats={'','','','0','','0.00000'};
    sheets=appendSheet(sheets,'ENZYMES',headers,enzSheet,false,enzColFormats);

    %ENZRXNS sheet: one row per enzyme-constrained reaction. The ENZYMES
    %column encodes the subunit stoichiometry from model.ec.rxnEnzMat as
    %'enzyme:count' pairs (e.g. 'P12345:1;P67890:2').
    headers={'#';'ID';'KCAT';'SOURCE';'NOTE';'EC-NUMBER';'ENZYMES'};
    ecRxnSheet=cell(numel(model.ec.rxns),numel(headers));
    ecRxnSheet(:,2)=model.ec.rxns(:);
    ecRxnSheet(:,3)=numToCell(model.ec.kcat);
    ecRxnSheet(:,4)=model.ec.source(:);
    ecRxnSheet(:,5)=model.ec.notes(:);
    ecRxnSheet(:,6)=model.ec.eccodes(:);
    ecRxnSheet(:,7)=ecEnzymePairs(model.ec);
    sheets=appendSheet(sheets,'ENZRXNS',headers,ecRxnSheet,false,{});
end

%Write all accumulated sheets to the .xlsx file in one pass.
writeExcel(fileName,sheets);
end

function sheets=appendSheet(sheets,name,header,data,isIntegers,colFormats)
%Append one worksheet definition for writeExcel (see sheetWidths for the
%per-sheet column widths).
n=numel(sheets)+1;
sheets(n).name=name;
sheets(n).header=header;
sheets(n).data=data;
sheets(n).widths=sheetWidths(name);
sheets(n).isIntegers=isIntegers;
sheets(n).colFormats=colFormats;
end

function w=sheetWidths(name)
%Per-sheet column widths. The stored values are in 1/256th of a character
%and are converted to character units below.
switch name
    case 'RXNS';    w256=[786;2358;7860;15719;3406;7860;3406;3406;2358;3406;7860;7860;3668;7860;7860;4192];
    case 'METS';    w256=[786;7860;7860;3668;7860;7860;7860;3406;3668;1834];
    case 'COMPS';   w256=[786;3144;7860;3144;7860];
    case 'GENES';   w256=[786;3144;7860;3144;3406];
    case 'MODEL';   w256=[786;3144;7860;3668;3668;5240;5240;7860;7860;2620;7860];
    case 'ENZYMES'; w256=[786;3144;3144;3406;15719;3406];
    case 'ENZRXNS'; w256=[786;5240;3406;5240;7860;5240;7860];
    otherwise;      w256=[];
end
w=w256(:).'/256;
end

function c=numToCell(v)
%Convert a numeric vector to a column cell array, representing NaN values
%as empty (blank) cells so they are written as empty cells in the sheet.
c=num2cell(v(:));
c(cellfun(@(x) isnumeric(x) && isnan(x),c))={[]};
end

function c=ecEnzymePairs(ec)
%Build the ENZRXNS 'ENZYMES' column. For each ec reaction, list the enzymes
%and their subunit stoichiometry from rxnEnzMat as 'enzyme:count;...'.
%Reactions without enzymes get a blank cell.
nRxns=numel(ec.rxns);
c=cell(nRxns,1);
if ~isfield(ec,'rxnEnzMat') || isempty(ec.rxnEnzMat)
    return
end
for i=1:min(nRxns,size(ec.rxnEnzMat,1))
    j=find(ec.rxnEnzMat(i,:));
    if isempty(j)
        continue
    end
    pairs=cell(1,numel(j));
    for k=1:numel(j)
        pairs{k}=[ec.enzymes{j(k)} ':' num2str(full(ec.rxnEnzMat(i,j(k))))];
    end
    c{i}=strjoin(pairs,';');
end
end
