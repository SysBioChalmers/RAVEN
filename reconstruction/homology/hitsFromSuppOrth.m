function [blastStructure,suppOrthInfo]=hitsFromSuppOrth(source,organismID,varargin)
% hitsFromSuppOrth  Read SuppOrth consensus orthologs into a blastStructure.
%
% SuppOrth (https://github.com/fma91/SuppOrth) runs several ortholog
% predictors (OrthoFinder, SonicParanoid2, Broccoli, and optionally an
% imported OMA run) on a set of proteomes and writes suppOrth.tsv: one row
% per cross-species protein pair, with the number of predictors that
% support it. This function turns that table into a blastStructure for
% getModelFromHomology.
%
% Parameters
% ----------
% source : char or table
%     path to suppOrth.tsv, or a MATLAB table with its columns. Protein
%     ids are read as text, so ids such as "NA" or "0123" are kept.
% organismID : char
%     SuppOrth species label of the organism to build a model for: the
%     FASTA file name without extension, or the --labels value. It becomes
%     fromId/toId unless speciesToModel maps it to another id.
%
% Name-Value Arguments
% --------------------
% speciesToModel : cell or containers.Map
%     SuppOrth species label to model id, as an N-by-2 cell array or a
%     containers.Map. A label mapped to '' is left out. When given, every
%     label in the table must be included. When empty, labels are used as
%     model ids (default []).
% minSupport : double
%     minimum number of predictors that found the pair (Lv_support).
%     SuppOrth has already filtered the table with its own --supportLv and
%     --identity_filter; this filter is applied in addition (default 1).
% requireTools : cell
%     predictors that must all have found the pair, by column name
%     ('orthofinder', 'sonicparanoid', 'broccoli', 'oma') (default {}).
% minIdentity : double
%     minimum SuppOrth identity. A pair without an alignment hit fails this
%     filter (default [], no filter).
% maxEvalue : double
%     maximum SuppOrth e-value. A pair without an alignment hit fails this
%     filter (default [], no filter).
% proteinToGene : cell, containers.Map or char
%     protein id to model gene id, as an N-by-2 cell array, a
%     containers.Map, a JSON file with one object of strings, or a
%     two-column tab-separated file (lines starting with # are skipped).
%     Proteins missing from the map keep their id, and their number is
%     reported in a warning. A gene pair takes the highest support over its
%     protein pairs, the union of their predictors, the highest identity
%     and bitscore, and the lowest e-value; these values can come from
%     different protein pairs (default [], no mapping).
%
% Returns
% -------
% blastStructure : struct
%     one element per search direction between organismID and each other
%     species, as from getBlast or getDiamond. Every pair is written in
%     both directions, so strictness 1 and 2 of getModelFromHomology give
%     the same result. evalue, identity, aligLen and ppos are placeholders
%     (0, 100, 1000, 100) that pass the filters of getModelFromHomology;
%     bitscore is SuppOrth's value (0 when it has none), which strictness 3
%     ranks by.
% suppOrthInfo : table
%     one row per gene pair, with modelId, newGene, templateGene, support,
%     supportedBy, identity, evalue and bitscore (SuppOrth's values).
%
% Notes
% -----
% Only pairs between organismID and another species are kept. Species that
% are not templates in getModelFromHomology have to be left out with
% speciesToModel, because getModelFromHomology requires every id in the
% blastStructure to be the new organism or a template.
%
% Examples
% --------
%     blastStructure = hitsFromSuppOrth('out/suppOrth.tsv','klac',...
%        'minSupport',2,'proteinToGene','klac_genes.tsv');
%
% See also
% --------
% runSuppOrth, getModelFromHomology, getDiamond

p=parseRAVENargs(varargin, {'speciesToModel',[]; 'minSupport',1; ...
    'requireTools',{}; 'minIdentity',[]; 'maxEvalue',[]; 'proteinToGene',[]});
organismID=char(organismID);
requireTools=cellstr(p.requireTools);
if isempty(p.requireTools)
    requireTools={};
end

T=readTable(source);
columns=T.Properties.VariableNames;

labels=unique([T.Species1;T.Species2]);
if ~ismember(organismID,labels)
    error('RAVEN:badInput','organismID "%s" is not a species in the table: %s',...
        organismID,strjoin(labels,', '));
end
if isempty(p.speciesToModel)
    idOf=containers.Map(labels,labels);
else
    idOf=toMap(p.speciesToModel);
    unmapped=labels(~isKey(idOf,labels));
    if ~isempty(unmapped)
        error('RAVEN:badInput',['speciesToModel has no entry for %s; map each '...
            'label to a model id, or to '''' to leave it out.'],strjoin(unmapped,', '));
    end
    if isempty(idOf(organismID))
        error('RAVEN:badInput','speciesToModel maps organismID "%s" to ''''.',organismID);
    end
end

%The per-predictor 0/1 columns are those between SupportedBy and Identity
alignment={'Identity','e_val','Bitscore'};
start=find(strcmp(columns,'SupportedBy'))+1;
stop=find(strcmp(columns,'Identity'))-1;
if isempty(stop)
    stop=numel(columns);
end
toolColumns=setdiff(columns(start:stop),alignment,'stable');
unknown=setdiff(requireTools,toolColumns,'stable');
if ~isempty(unknown)
    error('RAVEN:badInput','requireTools %s not in the table; its predictors are %s',...
        strjoin(unknown,', '),strjoin(toolColumns,', '));
end

keep=T.Lv_support>=p.minSupport;
for i=1:numel(requireTools)
    keep=keep & toNumber(T.(requireTools{i}))==1;
end
if ~isempty(p.minIdentity)
    keep=keep & T.Identity>=p.minIdentity; %NaN fails
end
if ~isempty(p.maxEvalue)
    keep=keep & T.e_val<=p.maxEvalue; %NaN fails
end
T=T(keep,:);

%Orient every pair as (organism protein, template species, template protein)
first=strcmp(T.Species1,organismID) & ~strcmp(T.Species2,organismID);
second=strcmp(T.Species2,organismID) & ~strcmp(T.Species1,organismID);
newProt=[T.Protein1(first);T.Protein2(second)];
label=[T.Species2(first);T.Species1(second)];
template=[T.Protein2(first);T.Protein1(second)];
support=[T.Lv_support(first);T.Lv_support(second)];
supportedBy=[T.SupportedBy(first);T.SupportedBy(second)];
identity=[T.Identity(first);T.Identity(second)];
evalue=[T.e_val(first);T.e_val(second)];
bitscore=[T.Bitscore(first);T.Bitscore(second)];
modelId=values(idOf,label);
modelId=modelId(:);
I=~cellfun(@isempty,modelId);
[newProt,template,modelId,support,supportedBy,identity,evalue,bitscore]=...
    deal(newProt(I),template(I),modelId(I),support(I),supportedBy(I),...
    identity(I),evalue(I),bitscore(I));

if ~isempty(p.proteinToGene)
    geneOf=toMap(p.proteinToGene);
    proteins=unique([newProt;template]);
    nMissing=sum(~isKey(geneOf,proteins));
    if nMissing>0
        warning('RAVEN:hitsFromSuppOrth:unmappedProteins',['%d of %d proteins are not '...
            'in proteinToGene and keep their protein id.'],nMissing,numel(proteins));
    end
    newProt=mapIds(geneOf,newProt);
    template=mapIds(geneOf,template);
end

if isempty(newProt)
    error('RAVEN:badInput',['No SuppOrth pairs between "%s" and a template species '...
        'remain after filtering.'],organismID);
end

%One row per gene pair, sorted by model, new gene and template gene
key=strcat(modelId,char(31),newProt,char(31),template);
[key,firstRow,group]=unique(key);
info=table();
info.modelId=modelId(firstRow);
info.newGene=newProt(firstRow);
info.templateGene=template(firstRow);
info.support=accumarray(group,support,[],@max);
rows=accumarray(group,(1:numel(group))',[],@(r) {r});
info.supportedBy=cellfun(@(r) joinLabels(supportedBy(r)),rows,'UniformOutput',false);
info.identity=accumarray(group,identity,[],@max);
info.evalue=accumarray(group,evalue,[],@min);
info.bitscore=accumarray(group,bitscore,[],@max);
suppOrthInfo=info;

newId=idOf(organismID);
score=info.bitscore;
score(isnan(score))=0;
models=unique(info.modelId);
blastStructure=[];
for direction=1:2
    for i=1:numel(models)
        r=strcmp(info.modelId,models{i});
        element=[];
        if direction==1
            element.fromId=models{i};
            element.toId=newId;
            element.fromGenes=info.templateGene(r);
            element.toGenes=info.newGene(r);
        else
            element.fromId=newId;
            element.toId=models{i};
            element.fromGenes=info.newGene(r);
            element.toGenes=info.templateGene(r);
        end
        n=sum(r);
        element.evalue=zeros(n,1);
        element.identity=100*ones(n,1);
        element.aligLen=1000*ones(n,1);
        element.bitscore=score(r);
        element.ppos=100*ones(n,1);
        blastStructure=[blastStructure element]; %#ok<AGROW>
    end
end
end

function T=readTable(source)
required={'Species1','Species2','Protein1','Protein2','Lv_support','SupportedBy'};
text=[required(1:4),{'SupportedBy'}];
if istable(source)
    T=source;
else
    source=char(source);
    if ~isfile(source)
        error('RAVEN:badInput','SuppOrth table "%s" does not exist.',source);
    end
    lines=splitlines(fileread(source));
    lines=lines(~cellfun(@isempty,lines));
    header=strsplit(lines{1},'\t','CollapseDelimiters',false);
    fields=regexp(lines(2:end),'\t','split');
    width=cellfun(@numel,fields);
    bad=find(width~=numel(header),1);
    if ~isempty(bad)
        error('RAVEN:badInput','%s: line %d has %d fields, the header has %d.',...
            source,bad+1,width(bad),numel(header));
    end
    fields=vertcat(fields{:});
    if isempty(fields)
        fields=cell(0,numel(header));
    end
    T=cell2table(fields,'VariableNames',header);
end
missing=setdiff(required,T.Properties.VariableNames,'stable');
if ~isempty(missing)
    error('RAVEN:badInput',['Not a suppOrth.tsv table (SuppOrth 0.4): missing columns %s. '...
        'suppOrth.genes.tsv is not supported; use proteinToGene instead.'],strjoin(missing,', '));
end
%Protein ids such as "NA" or "0123" stay text
for i=1:numel(text)
    T.(text{i})=cellstr(string(T.(text{i})));
end
for c={'Identity','e_val','Bitscore'}
    if ismember(c{1},T.Properties.VariableNames)
        T.(c{1})=toNumber(T.(c{1}));
    else
        T.(c{1})=nan(height(T),1);
    end
end
T.Lv_support=toNumber(T.Lv_support);
end

function x=toNumber(x)
%Read numbers exactly; empty and non-numeric fields become NaN
if ~isnumeric(x)
    x=str2double(x);
end
x=double(x(:));
end

function m=toMap(x)
if isa(x,'containers.Map')
    m=x;
elseif iscell(x)
    if size(x,2)~=2
        error('RAVEN:badInput','A cell array map needs two columns.');
    end
    m=containers.Map(cellstr(string(x(:,1))),cellstr(string(x(:,2))));
elseif ischar(x) || isstring(x)
    m=readMapFile(char(x));
else
    error('RAVEN:badInput','A map must be a two-column cell array, a containers.Map or a file.');
end
end

function m=readMapFile(file)
[~,~,ext]=fileparts(file);
if any(strcmpi(ext,{'.pkl','.pickle'}))
    error('RAVEN:badInput',['%s: pickle maps are not read; use a cell array, '...
        'a JSON file or a two-column tab-separated file.'],file);
end
if ~isfile(file)
    error('RAVEN:badInput','Map file "%s" does not exist.',file);
end
text=fileread(file);
if strcmpi(ext,'.json')
    %Keys are protein ids, which jsondecode would turn into field names
    pairs=regexp(text,'"((?:[^"\\]|\\.)*)"\s*:\s*"((?:[^"\\]|\\.)*)"','tokens');
    pairs=vertcat(pairs{:});
    if isempty(pairs)
        m=containers.Map('KeyType','char','ValueType','char');
        return
    end
    pairs=regexprep(pairs,'\\(["\\/])','$1');
    m=containers.Map(pairs(:,1),pairs(:,2));
    return
end
lines=splitlines(text);
lines=lines(~cellfun(@isempty,strtrim(lines)) & ~startsWith(lines,'#'));
fields=regexp(lines,'\t','split');
short=find(cellfun(@numel,fields)<2,1);
if ~isempty(short)
    error('RAVEN:badInput','%s: expected two tab-separated columns, got "%s"',file,lines{short});
end
keys=strtrim(cellfun(@(f) f{1},fields,'UniformOutput',false));
vals=strtrim(cellfun(@(f) f{2},fields,'UniformOutput',false));
if isempty(keys)
    m=containers.Map('KeyType','char','ValueType','char');
else
    m=containers.Map(keys,vals);
end
end

function ids=mapIds(m,ids)
known=isKey(m,ids);
if any(known)
    ids(known)=values(m,ids(known));
end
end

function joined=joinLabels(labels)
%Union of comma-separated predictor labels, in SuppOrth's F, S, B, O order
found=strtrim(split(strjoin(labels,','),','));
found=unique(found(~cellfun(@isempty,found)));
order={'F','S','B','O'};
[isKnown,rank]=ismember(found,order);
known=found(isKnown);
[~,o]=sort(rank(isKnown));
joined=strjoin([known(o);sort(found(~isKnown))],',');
end
