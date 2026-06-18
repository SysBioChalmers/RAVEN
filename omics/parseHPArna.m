function arrayData=parseHPArna(fileName, varargin)
% parseHPArna  Parse a dump of Human Protein Atlas (HPA) RNA-Seq data.
%
% Parameters
% ----------
% fileName : char
%     tab-separated database dump of HPA RNA data. For details regarding the
%     format, see http://www.proteinatlas.org/about/download.
%
% Name-Value Arguments
% --------------------
% version : double
%     accepted for backward compatibility but ignored; the format is now
%     inferred from the column headers (default 19).
%
% Returns
% -------
% arrayData : struct
%     parsed HPA RNA data with fields:
%
%     - genes : cell array with the unique Ensembl gene IDs
%     - geneNames : cell array with the gene symbols
%     - tissues : cell array with the unique tissue (or sample) names
%     - levels : matrix of expression values (rows = genes, cols = tissues)
%
% Examples
% --------
%     arrayData = parseHPArna(fileName);

p=parseRAVENargs(varargin, {'version',19}); %#ok<NASGU> retained for backward compat

fileName=char(fileName);
if ~isfile(fileName)
    error('HPA file %s cannot be found', string(fileName));
end

%Parse header and all data columns (RNA files are always TSV)
fid=fopen(fileName,'r');
headers=strsplit(fgetl(fid), sprintf('\t'));
fmtStr=strip(repmat('%q ',1,numel(headers)));
hpa=textscan(fid, fmtStr, 'Delimiter','\t');
fclose(fid);

%Helper: column index by name (returns empty if not found)
ci=@(name) find(strcmp(headers,name),1);

%Validate required columns
for h={'Gene','Gene name'}
    if isempty(ci(h{1}))
        error(['Could not find the column "' h{1} '". ' ...
            'Make sure the input file matches the format at https://www.proteinatlas.org/about/download']);
    end
end

%Detect tissue column: 'Tissue' (v19) or 'Sample' (v18)
tissueCol='Tissue';
if isempty(ci('Tissue')) && ~isempty(ci('Sample'))
    tissueCol='Sample';
elseif isempty(ci('Tissue'))
    error(['Could not find a tissue column ("Tissue" or "Sample"). ' ...
        'Make sure the input file matches the format at https://www.proteinatlas.org/about/download']);
end

%Detect expression value column: 'TPM' (v19) or 'Value' (v18)
valueCol='TPM';
if isempty(ci('TPM')) && ~isempty(ci('Value'))
    valueCol='Value';
elseif isempty(ci('TPM'))
    error(['Could not find an expression value column ("TPM" or "Value"). ' ...
        'Make sure the input file matches the format at https://www.proteinatlas.org/about/download']);
end

[arrayData.genes, P, I]=unique(hpa{ci('Gene')});
arrayData.geneNames=hpa{ci('Gene name')}(P);
[arrayData.tissues, ~, J]=unique(hpa{ci(tissueCol)});

arrayData.levels=NaN(max(I),max(J));
linearInd=sub2ind(size(arrayData.levels),I,J);
arrayData.levels(linearInd)=str2double(hpa{ci(valueCol)});
end
