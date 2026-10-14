function [blastStructure,suppOrthInfo]=runSuppOrth(organismID,fastaFile,modelIDs,refFastaFiles,outDir,varargin)
% runSuppOrth  Run SuppOrth on an organism and template proteomes.
%
% Runs the suppOrth command (https://github.com/fma91/SuppOrth) and reads
% its consensus orthologs with hitsFromSuppOrth.
%
% Parameters
% ----------
% organismID : char
%     the id of the organism of interest, as in getModelFromHomology.
% fastaFile : char
%     a FASTA file (optionally gzipped) with the proteins of the organism.
% modelIDs : cell
%     the template model ids, as "model.id" in getModelFromHomology.
% refFastaFiles : cell
%     the proteome of each template, in the order of modelIDs, with FASTA
%     ids equal to the gene ids of that model (or mapped to them with
%     proteinToGene).
% outDir : char
%     folder for SuppOrth's output.
%
% Name-Value Arguments
% --------------------
% tools : char
%     comma-separated predictors, or 'all' for orthofinder, sonicparanoid
%     and broccoli (default 'all').
% supportLv : double
%     SuppOrth's --supportLv: a pair is kept if at least this many
%     predictors found it, or if its identity is above identityFilter
%     (default 3, or the number of tools when fewer are selected).
% identityFilter : double
%     SuppOrth's --identity_filter (default 80).
% alignment : char
%     SuppOrth's --alignment (default 'diamond').
% threads : double
%     number of threads (default: the cores assigned to MATLAB).
% extraArgs : cell
%     further arguments for "suppOrth run", for example {'--resume'} or
%     {'--oma','path/to/Output'} (default {}).
% suppOrth : char
%     the suppOrth command. If empty, the RAVEN_SUPPORTH environment
%     variable is used, and then suppOrth on the PATH (default '').
%
% The name-value arguments of hitsFromSuppOrth (minSupport, requireTools,
% minIdentity, maxEvalue, proteinToGene, speciesToModel) are passed on.
%
% Returns
% -------
% blastStructure : struct
%     as from hitsFromSuppOrth, for getModelFromHomology.
% suppOrthInfo : table
%     SuppOrth's values per gene pair, as from hitsFromSuppOrth.
%
% Notes
% -----
% SuppOrth requires protein ids to be unique over all proteomes. When two
% FASTA files share an id (for example two templates of the same species),
% every id is given a "<label>__" prefix in copies written to
% outDir/prefixed_inputs, and the prefix is removed from the returned
% hits. outDir/suppOrth.tsv then contains the prefixed ids. SuppOrth's
% console output is written to outDir/run_supporth.log.
%
% SuppOrth's predictors run on Linux and macOS only. On Windows, run
% SuppOrth under WSL or on a Linux host, and read suppOrth.tsv with
% hitsFromSuppOrth.
%
% Examples
% --------
%     blastStructure = runSuppOrth('klac','klac.fasta',{'yeastGEM'},...
%        {'yeast.fasta'},'supporth_out','proteinToGene','klac_genes.tsv');
%     draftModel = getModelFromHomology(yeastGEM,blastStructure,'klac');
%
% See also
% --------
% hitsFromSuppOrth, getModelFromHomology, getDiamond

readerArgs={'speciesToModel',[]; 'minSupport',1; 'requireTools',{}; ...
    'minIdentity',[]; 'maxEvalue',[]; 'proteinToGene',[]};
p=parseRAVENargs(varargin, [{'tools','all'; 'supportLv',[]; 'identityFilter',80; ...
    'alignment','diamond'; 'threads',[]; 'extraArgs',{}; 'suppOrth',''}; readerArgs]);

if ispc
    error('RAVEN:badInput',['SuppOrth''s predictors run on Linux and macOS only. Run '...
        'SuppOrth under WSL or on a Linux host, and read suppOrth.tsv with hitsFromSuppOrth.']);
end
organismID=char(organismID);
modelIDs=cellstr(modelIDs);
refFastaFiles=cellstr(refFastaFiles);
if numel(modelIDs)~=numel(refFastaFiles)
    error('RAVEN:badInput','modelIDs and refFastaFiles must have the same length.');
end
if isempty(modelIDs)
    error('RAVEN:badInput','At least one template is needed.');
end
labels=[{organismID};modelIDs(:)];
for i=1:numel(labels)
    l=labels{i};
    if isempty(l) || any(strcmp(l,{'.','..'})) || ~strcmp(l,strtrim(l)) || any(ismember(l,',/\'))
        error('RAVEN:badInput',['"%s" cannot be a SuppOrth species label: it must be '...
            'non-empty, not "." or "..", have no leading or trailing space, and contain '...
            'no ",", "/" or "\\".'],l);
    end
end
if numel(unique(labels))~=numel(labels)
    error('RAVEN:badInput','organismID and modelIDs must be distinct, got %s',strjoin(labels,', '));
end

tools=strtrim(strsplit(char(p.tools),','));
tools=tools(~cellfun(@isempty,tools));
if isequal(tools,{'all'})
    tools={'orthofinder','sonicparanoid','broccoli'};
end
supportLv=p.supportLv;
if isempty(supportLv)
    supportLv=min(3,numel(tools));
end
if supportLv<1 || supportLv>4
    error('RAVEN:badInput','supportLv must be between 1 and 4, got %d',supportLv);
end
extraArgs=cellstr(p.extraArgs);
if isempty(p.extraArgs)
    extraArgs={};
end
if supportLv>numel(tools) && ~any(strcmp(extraArgs,'--oma'))
    warning('RAVEN:runSuppOrth:supportLv',['supportLv=%d is more than the %d predictor(s) '...
        'selected; SuppOrth then keeps only pairs above identityFilter (%g %%).'],...
        supportLv,numel(tools),p.identityFilter);
end
threads=p.threads;
if isempty(threads)
    cores=evalc('feature(''numcores'')');
    cores=strsplit(cores,'MATLAB was assigned: ');
    threads=str2double(regexp(cores{2},'^\d*','match','once'));
end

exe=resolveSuppOrth(char(p.suppOrth));
fastas=[{char(fastaFile)};refFastaFiles(:)];
for i=1:numel(fastas)
    if ~isfile(fastas{i})
        error('RAVEN:badInput','FASTA file "%s" does not exist.',fastas{i});
    end
end
if ~isfolder(outDir)
    mkdir(outDir);
end

%SuppOrth requires protein ids to be unique over all proteomes. Templates of
%the same species often share gene ids, so in that case every id gets a
%"<label>__" prefix in a copy of its FASTA, and the prefix is removed below.
ids=cell(numel(fastas),1);
for i=1:numel(fastas)
    ids{i}=fastaIds(fastas{i});
    ids{i}=unique(ids{i});
end
allIds=vertcat(ids{:});
shared=numel(unique(allIds))<numel(allIds);
if shared
    inputs=fullfile(outDir,'prefixed_inputs');
    if ~isfolder(inputs)
        mkdir(inputs);
    end
    for i=1:numel(fastas)
        dest=fullfile(inputs,[labels{i} '.faa']);
        writePrefixed(fastas{i},dest,[labels{i} '__']);
        fastas{i}=dest;
    end
end

args=[{exe,'run'},fastas',{'--labels',strjoin(labels,','),'-o',char(outDir),...
    '--tools',strjoin(tools,','),'--supportLv',num2str(supportLv),...
    '--identity_filter',num2str(p.identityFilter),'--alignment',char(p.alignment),...
    '-t',num2str(threads)},extraArgs(:)'];
cmd=strjoin(cellfun(@shellQuote,args,'UniformOutput',false),' ');
%MATLAB's library path can break the predictors' compiled dependencies
[status,output]=system(['env -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH ' cmd ' 2>&1']);
logFile=fullfile(outDir,'run_supporth.log');
fid=fopen(logFile,'w');
fprintf(fid,'%s\n\n%s\n',strjoin(args,' '),strtrim(output));
fclose(fid);
if status~=0
    lines=splitlines(strtrim(output));
    %SuppOrth reports a failed stage before the consensus summary, which can
    %push it out of the last lines
    head=lines(1:max(0,end-20));
    failed=head(contains(head,'FAILED') | startsWith(head,'error'));
    tail=[failed;lines(max(1,end-19):end)];
    error('RAVEN:runSuppOrth:failed',['suppOrth exited with code %d; its output in %s '...
        'is incomplete and was not read. Full output in %s.\n%s'],status,char(outDir),...
        logFile,strjoin(tail,newline));
end

source=fullfile(outDir,'suppOrth.tsv');
if shared
    source=stripPrefixes(source);
end
readerValues=cellfun(@(n) p.(n),readerArgs(:,1),'UniformOutput',false);
reader=[readerArgs(:,1),readerValues]';
[blastStructure,suppOrthInfo]=hitsFromSuppOrth(source,organismID,reader{:});
end

function exe=resolveSuppOrth(exe)
if ~isempty(exe)
    return
end
exe=getenv('RAVEN_SUPPORTH');
if ~isempty(exe)
    return
end
for name={'suppOrth','supporth'}
    [status,found]=system(['command -v ' name{1}]);
    if status==0 && ~isempty(strtrim(found))
        exe=strtrim(found);
        return
    end
end
error('RAVEN:runSuppOrth:notFound',['Could not find the SuppOrth command suppOrth. '...
    'Install it with pip install ''git+https://github.com/fma91/SuppOrth.git@stand_alone'' followed '...
    'by suppOrth install, then put suppOrth on the PATH, set RAVEN_SUPPORTH, or pass '...
    'the suppOrth argument.']);
end

function lines=readFastaLines(file)
[~,~,ext]=fileparts(file);
if strcmpi(ext,'.gz')
    tmp=tempname;
    unzipped=gunzip(file,tmp);
    text=fileread(unzipped{1});
    rmdir(tmp,'s');
else
    text=fileread(file);
end
lines=splitlines(text);
end

function ids=fastaIds(file)
%Sequence ids as SuppOrth reads them: the first word of each header line
lines=readFastaLines(file);
headers=lines(startsWith(lines,'>'));
words=regexp(headers,'^>\s*(\S+)','tokens','once');
words=words(~cellfun(@isempty,words));
ids=cellfun(@(w) w{1},words,'UniformOutput',false);
end

function writePrefixed(src,dest,prefix)
lines=readFastaLines(src);
isHeader=startsWith(lines,'>');
lines(isHeader)=strcat('>',prefix,extractAfter(lines(isHeader),1));
if ~isempty(lines) && isempty(lines{end})
    lines(end)=[];
end
fid=fopen(dest,'w');
fprintf(fid,'%s\n',lines{:});
fclose(fid);
end

function T=stripPrefixes(file)
%Remove the "<label>__" prefix that runSuppOrth added to protein ids
lines=splitlines(fileread(file));
lines=lines(~cellfun(@isempty,lines));
header=strsplit(lines{1},'\t','CollapseDelimiters',false);
fields=regexp(lines(2:end),'\t','split');
fields=vertcat(fields{:});
if isempty(fields)
    fields=cell(0,numel(header));
end
T=cell2table(fields,'VariableNames',header);
for pair={{'Species1','Protein1'},{'Species2','Protein2'}}
    s=T.(pair{1}{1});
    prot=T.(pair{1}{2});
    prefix=strcat(s,'__');
    has=cellfun(@startsWith,prot,prefix);
    prot(has)=cellfun(@(x,pre) x(numel(pre)+1:end),prot(has),prefix(has),'UniformOutput',false);
    T.(pair{1}{2})=prot;
end
end

function q=shellQuote(s)
q=['''' strrep(char(s),'''','''\''''') ''''];
end
