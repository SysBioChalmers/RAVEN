function hpaData=parseHPA(fileName, varargin)
% parseHPA  Parse a database dump of the Human Protein Atlas (HPA).
%
% Parameters
% ----------
% fileName : char
%     comma- or tab-separated database dump of HPA protein data. For details
%     regarding the format, see http://www.proteinatlas.org/about/download.
%
% Name-Value Arguments
% --------------------
% version : double
%     accepted for backward compatibility but ignored; the format is now
%     inferred from the column headers (default 19).
%
% Returns
% -------
% hpaData : struct
%     parsed HPA data with fields:
%
%     - genes : cell array with the unique gene IDs (Ensembl in v>=18)
%     - geneNames : cell array with the gene symbols (present in v>=18 only)
%     - tissues : cell array with the tissue names (may not be unique; one
%       entry per tissue-cell-type combination)
%     - celltypes : cell array with the cell type names for each tissue
%     - levels : cell array with the unique expression level labels
%     - types : cell array with the unique evidence types (v<18 only)
%     - reliabilities : cell array with the unique reliability levels
%     - gene2Level : sparse gene × tissue-cell-type matrix; value i,j is
%       the index in hpaData.levels for gene i in tissue-cell-type j
%     - gene2Type : sparse gene × tissue-cell-type matrix (v<18 only)
%     - gene2Reliability : sparse gene × tissue-cell-type matrix; value i,j
%       is the index in hpaData.reliabilities for gene i in tissue-cell-type j
%
% Examples
% --------
%     hpaData = parseHPA(fileName);

p=parseRAVENargs(varargin, {'version',19}); %#ok<NASGU> retained for backward compat

fileName=char(fileName);
if ~isfile(fileName)
    error('HPA file %s cannot be found', string(fileName));
end

%Read header line and auto-detect delimiter
fid=fopen(fileName,'r');
headerLine=fgetl(fid);
fclose(fid);

if any(headerLine==9)  % char(9) = horizontal tab
    delim=sprintf('\t');
else
    delim=',';
end

%Parse header and all data columns
fid=fopen(fileName,'r');
headers=strsplit(fgetl(fid), delim);
fmtStr=strip(repmat('%q ',1,numel(headers)));
hpa=textscan(fid, fmtStr, 'Delimiter',delim);
fclose(fid);

%Helper: column index by name
ci=@(name) find(strcmp(headers,name),1);

%Validate required columns common to all versions
for h={'Gene','Tissue','Cell type','Level','Reliability'}
    if isempty(ci(h{1}))
        error(['Could not find the column "' h{1} '". ' ...
            'Make sure the input file matches the format at https://www.proteinatlas.org/about/download']);
    end
end

hasGeneName=~isempty(ci('Gene name'));
hasExprType=~isempty(ci('Expression type'));

[hpaData.genes, P, I]=unique(hpa{ci('Gene')});
if hasGeneName
    hpaData.geneNames=hpa{ci('Gene name')}(P);
end

tissueCellStr=strcat(hpa{ci('Tissue')},'€',hpa{ci('Cell type')});
[~, J, K]=unique(tissueCellStr);
hpaData.tissues=hpa{ci('Tissue')}(J);
hpaData.celltypes=hpa{ci('Cell type')}(J);

[hpaData.levels, ~, L]=unique(hpa{ci('Level')});
[hpaData.reliabilities, ~, N]=unique(hpa{ci('Reliability')});

nGenes=numel(hpaData.genes);
nTissues=max(K);

hpaData.gene2Level=sparse(I,K,L,nGenes,nTissues);
hpaData.gene2Reliability=sparse(I,K,N,nGenes,nTissues);

if hasExprType
    [hpaData.types, ~, M]=unique(hpa{ci('Expression type')});
    hpaData.gene2Type=sparse(I,K,M,nGenes,nTissues);
end
end
