function writeExcel(fileName,sheets)
% writeExcel  Write one or more sheets to an .xlsx (Office Open XML) file.
%
% Writes spreadsheet data to the SpreadsheetML (.xlsx) format without any
% external library (no Apache POI / Java) and without MATLAB toolboxes. An
% .xlsx file is a ZIP archive of XML parts; this function generates those
% parts and packs them with the built-in zip function.
%
% Parameters
% ----------
% fileName : char
%     path of the .xlsx file to write. An existing file is overwritten.
% sheets : struct
%     struct array, one element per worksheet, with fields:
%       name : char
%           the worksheet name (truncated to 31 characters; the characters
%           []*?/\: are removed).
%       data : cell
%           cell matrix with the cell values. Supported value types are char
%           (text), numeric scalar, logical, and empty (written as a blank
%           cell).
%       header : cell, optional
%           column captions written as a bold, frozen first row above the
%           data (default none).
%       widths : double, optional
%           vector of column widths in character units (default: column
%           widths are left to the spreadsheet application).
%       isIntegers : logical, optional
%           if true numeric cells use the '0' number format, otherwise
%           '##0.00' (default true).
%       colFormats : cell, optional
%           per-column number-format codes that override the sheet default
%           for individual columns (default none).
%
% Examples
% --------
%     s.name='SHEET1'; s.header={'ID','VALUE'}; s.data={'a',1;'b',2};
%     writeExcel('out.xlsx',s);
%
% See also
% --------
% exportToExcelFormat

%Normalise the sheet structs and collect the set of number-format codes used
%across all sheets, so each distinct format gets a single style entry.
formatCodes={};
S=struct('name',{},'header',{},'data',{},'widths',{},'defFmt',{},'colFmt',{});
for k=1:numel(sheets)
    in=sheets(k);
    S(k).name=cleanSheetName(in.name);
    if isfield(in,'header'); S(k).header=in.header; else; S(k).header={}; end
    S(k).data=in.data;
    if isfield(in,'widths'); S(k).widths=in.widths; else; S(k).widths=[]; end
    if isfield(in,'isIntegers') && ~isempty(in.isIntegers) && ~in.isIntegers
        S(k).defFmt='##0.00';
    else
        S(k).defFmt='0';
    end
    if isfield(in,'colFormats'); S(k).colFmt=in.colFormats; else; S(k).colFmt={}; end

    nCols=size(S(k).data,2);
    if ~isempty(S(k).header); nCols=max(nCols,numel(S(k).header)); end
    for c=1:nCols
        if c<=numel(S(k).colFmt) && ~isempty(S(k).colFmt{c})
            fc=S(k).colFmt{c};
        else
            fc=S(k).defFmt;
        end
        if ~any(strcmp(fc,formatCodes)); formatCodes{end+1}=fc; end %#ok<AGROW>
    end
end

%Build the parts in a temporary directory, then zip them into the .xlsx.
tmpDir=tempname;
mkdir(tmpDir);
mkdir(fullfile(tmpDir,'_rels'));
mkdir(fullfile(tmpDir,'xl'));
mkdir(fullfile(tmpDir,'xl','_rels'));
mkdir(fullfile(tmpDir,'xl','worksheets'));

writePart(fullfile(tmpDir,'[Content_Types].xml'),contentTypes(numel(S)));
writePart(fullfile(tmpDir,'_rels','.rels'),rootRels());
writePart(fullfile(tmpDir,'xl','workbook.xml'),workbookXml(S));
writePart(fullfile(tmpDir,'xl','_rels','workbook.xml.rels'),workbookRels(numel(S)));
writePart(fullfile(tmpDir,'xl','styles.xml'),stylesXml(formatCodes));
for k=1:numel(S)
    writePart(fullfile(tmpDir,'xl','worksheets',sprintf('sheet%d.xml',k)), ...
        worksheetXml(S(k),formatCodes));
end

%Pack the parts. zip() preserves the relative paths when called from inside
%the staging directory, which is exactly the layout an .xlsx requires.
zipTmp=[tempname '.zip'];
oldDir=cd(tmpDir);
cleanup=onCleanup(@()cd(oldDir));
zip(zipTmp,{'[Content_Types].xml','_rels','xl'});
clear cleanup
if isfile(fileName); delete(fileName); end
movefile(zipTmp,fileName);
rmdir(tmpDir,'s');
end

% -------------------------------------------------------------------------

function xml=worksheetXml(s,formatCodes)
data=s.data;
header=s.header;
hasHeader=~isempty(header);
nRows=size(data,1);
nCols=size(data,2);
if hasHeader; nCols=max(nCols,numel(header)); end

%Per-column style index into cellXfs (numeric/logical cells only); strings
%use the default style 0 and headers the bold style 1.
colStyle=zeros(1,nCols);
for c=1:nCols
    if c<=numel(s.colFmt) && ~isempty(s.colFmt{c}); fc=s.colFmt{c}; else; fc=s.defFmt; end
    colStyle(c)=1+find(strcmp(fc,formatCodes),1); %0=default,1=bold,2..=formats
end

rows=cell(1,hasHeader+nRows);
ri=0;
if hasHeader
    ri=ri+1;
    cells=cell(1,numel(header));
    for c=1:numel(header)
        cells{c}=sprintf('<c r="%s1" t="inlineStr" s="1"><is><t xml:space="preserve">%s</t></is></c>', ...
            colLetter(c),xmlEscape(toText(header{c})));
    end
    rows{ri}=['<row r="1">' strjoin(cells,'') '</row>'];
end
for r=1:nRows
    rnum=r+hasHeader;
    cells=cell(1,nCols);
    m=0;
    for c=1:nCols
        v=data{r,c};
        if iscell(v) && ~isempty(v); v=v{1}; end %tolerate {value,list} cells
        if isempty(v); continue; end
        ref=sprintf('%s%d',colLetter(c),rnum);
        if ischar(v) || isstring(v)
            cellStr=sprintf('<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>', ...
                ref,xmlEscape(char(v)));
        elseif islogical(v)
            cellStr=sprintf('<c r="%s" t="b" s="%d"><v>%d</v></c>',ref,colStyle(c),v(1));
        elseif isnumeric(v)
            if ~isscalar(v) || isnan(v); continue; end
            cellStr=sprintf('<c r="%s" s="%d"><v>%.15g</v></c>',ref,colStyle(c),v);
        else
            cellStr=sprintf('<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>', ...
                ref,xmlEscape(toText(v)));
        end
        m=m+1;
        cells{m}=cellStr;
    end
    rows{ri+r}=sprintf('<row r="%d">%s</row>',rnum,strjoin(cells(1:m),''));
end

lastRow=max(hasHeader+nRows,1);
dimRef=sprintf('A1:%s%d',colLetter(max(nCols,1)),lastRow);

if hasHeader
    views=['<sheetViews><sheetView workbookViewId="0">' ...
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>' ...
        '<selection pane="bottomLeft" activeCell="A2" sqref="A2"/>' ...
        '</sheetView></sheetViews>'];
else
    views='<sheetViews><sheetView workbookViewId="0"/></sheetViews>';
end

cols='';
if ~isempty(s.widths)
    parts=cell(1,numel(s.widths));
    for c=1:numel(s.widths)
        parts{c}=sprintf('<col min="%d" max="%d" width="%.2f" customWidth="1"/>',c,c,s.widths(c));
    end
    cols=['<cols>' strjoin(parts,'') '</cols>'];
end

xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"' ...
    ' xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">' ...
    sprintf('<dimension ref="%s"/>',dimRef) views cols ...
    '<sheetData>' strjoin(rows,'') '</sheetData></worksheet>'];
end

% -------------------------------------------------------------------------

function xml=stylesXml(formatCodes)
%Custom number formats start at id 164 (ids below are reserved/built-in).
numFmts='';
for i=1:numel(formatCodes)
    numFmts=[numFmts sprintf('<numFmt numFmtId="%d" formatCode="%s"/>', ...
        163+i,xmlEscape(formatCodes{i}))]; %#ok<AGROW>
end
if ~isempty(numFmts)
    numFmts=sprintf('<numFmts count="%d">%s</numFmts>',numel(formatCodes),numFmts);
end

%cellXfs: 0=default (general,normal), 1=bold header, 2..=number formats.
xfs='<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>';
xfs=[xfs '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'];
for i=1:numel(formatCodes)
    xfs=[xfs sprintf(['<xf numFmtId="%d" fontId="0" fillId="0" borderId="0"' ...
        ' xfId="0" applyNumberFormat="1"/>'],163+i)]; %#ok<AGROW>
end

xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">' ...
    numFmts ...
    '<fonts count="2">' ...
    '<font><sz val="11"/><name val="Calibri"/></font>' ...
    '<font><b/><sz val="11"/><name val="Calibri"/></font>' ...
    '</fonts>' ...
    '<fills count="2">' ...
    '<fill><patternFill patternType="none"/></fill>' ...
    '<fill><patternFill patternType="gray125"/></fill>' ...
    '</fills>' ...
    '<borders count="1"><border/></borders>' ...
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' ...
    sprintf('<cellXfs count="%d">',2+numel(formatCodes)) xfs '</cellXfs>' ...
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' ...
    '</styleSheet>'];
end

function xml=workbookXml(S)
tags=cell(1,numel(S));
for k=1:numel(S)
    tags{k}=sprintf('<sheet name="%s" sheetId="%d" r:id="rId%d"/>', ...
        xmlEscape(S(k).name),k,k);
end
xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"' ...
    ' xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">' ...
    '<sheets>' strjoin(tags,'') '</sheets></workbook>'];
end

function xml=workbookRels(nSheets)
rels=cell(1,nSheets+1);
for k=1:nSheets
    rels{k}=sprintf(['<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org' ...
        '/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>'],k,k);
end
rels{nSheets+1}=sprintf(['<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org' ...
    '/officeDocument/2006/relationships/styles" Target="styles.xml"/>'],nSheets+1);
xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' ...
    strjoin(rels,'') '</Relationships>'];
end

function xml=contentTypes(nSheets)
overrides=cell(1,nSheets);
for k=1:nSheets
    overrides{k}=sprintf(['<Override PartName="/xl/worksheets/sheet%d.xml"' ...
        ' ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'],k);
end
xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' ...
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' ...
    '<Default Extension="xml" ContentType="application/xml"/>' ...
    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' ...
    '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' ...
    strjoin(overrides,'') '</Types>'];
end

function xml=rootRels()
xml=['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' ...
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument"' ...
    ' Target="xl/workbook.xml"/></Relationships>'];
end

% -------------------------------------------------------------------------

function writePart(path,str)
fid=fopen(path,'w','n','UTF-8');
if fid==-1; error('writeExcel: could not open ''%s'' for writing.',path); end
fwrite(fid,str,'char');
fclose(fid);
end

function s=xmlEscape(s)
s=strrep(s,'&','&amp;');
s=strrep(s,'<','&lt;');
s=strrep(s,'>','&gt;');
%Strip characters that are illegal in XML 1.0 (keep tab, LF, CR).
s=regexprep(s,'[\x00-\x08\x0B\x0C\x0E-\x1F]','');
end

function t=toText(v)
if ischar(v)
    t=v;
elseif isstring(v)
    t=char(v);
elseif islogical(v)
    if v(1); t='TRUE'; else; t='FALSE'; end
else
    t=num2str(v);
end
end

function L=colLetter(n)
L='';
while n>0
    r=mod(n-1,26);
    L=[char('A'+r) L]; %#ok<AGROW>
    n=floor((n-1)/26);
end
end

function name=cleanSheetName(name)
name=regexprep(name,'[\[\]\*\?/\\:]','');
if numel(name)>31; name=name(1:31); end
end
