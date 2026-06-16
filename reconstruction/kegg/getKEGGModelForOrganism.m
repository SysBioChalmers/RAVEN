function model=getKEGGModelForOrganism(organismID,varargin)
% getKEGGModelForOrganism  Reconstruct a model from KEGG protein homology.
%
% Reconstructs a genome-scale metabolic model based on protein homology to
% the orthologies in KEGG. If the target species is not available in KEGG,
% the user must select a closely related species. It is also possible to
% circumvent protein homology search (see fastaFile parameter for more
% details).
%
% Parameters
% ----------
% organismID : char
%     three or four letter abbreviation of the organism (as used in KEGG).
%     If not available, use a closely related species. This is used for
%     determing the phylogenetic distance. Use 'eukaryotes' or 'prokaryotes'
%     to get a model for the whole domain. Only applicable if fastaFile is
%     empty, i.e. no homology search should be performed.
%
% Name-Value Arguments
% --------------------
% fastaFile : char
%     a FASTA file that contains the protein sequences of the organism for
%     which to reconstruct a model. If no FASTA file is supplied then a
%     model is reconstructed based only on the organism abbreviation. This
%     option ignores all settings except for keepSpontaneous,
%     keepUndefinedStoich, keepIncomplete and keepGeneral.
% dataDir : char
%     directory for which to retrieve the input data, styled as
%     kegg118_prokaryotes or kegg118_eukaryotes, indicating the KEGG version
%     and whether the HMMs were trained on pro- or eukaryotic sequences. The
%     directory name matches the published HMM library it is paired with. The
%     prebuilt concatenated KO HMM library (dataDir.hmm) is downloaded here
%     from the corresponding RAVEN release if not already present. May also
%     contain a dataDir\keggdb sub-folder with a local KEGG FTP dump, used to
%     build the global KEGG model. This parameter should ALWAYS be provided.
% outDir : char
%     directory to save the results from the quering of the hidden Markov
%     models. The output is specific for the input sequences and the
%     settings used. It is stored in this manner so that the function can
%     continue if interrupted or if it should run in parallel. Be careful
%     not to leave output files from different organisms or runs with
%     different settings in the same folder. They will not be overwritten
%     (default is a temporary dir where all *.out files are deleted before
%     and after doing the reconstruction).
% keepSpontaneous : logical
%     include reactions labeled as "spontaneous" (default true).
% keepUndefinedStoich : logical
%     include reactions in the form n A <=> n+1 A. These will be dealt with
%     as two separate metabolites (default true).
% keepIncomplete : logical
%     include reactions which have been labelled as "incomplete",
%     "erroneous" or "unclear" (default true).
% keepGeneral : logical
%     include reactions which have been labelled as "general reaction".
%     These are reactions on the form "an aldehyde <=> an alcohol", and are
%     therefore unsuited for modelling purposes. Note that not all reactions
%     have this type of annotation, and the script will therefore not be
%     able to remove all such reactions (default false).
% cutOff : double
%     significance score from HMMer needed to assign genes to a KO (default
%     10^-50).
% minScoreRatioKO : double
%     ignore genes in a KO if their score is <log(score)/log(best score in
%     KO). This is to "prune" KOs which have many genes and where some are
%     clearly a better fit (default 0.3, lower is less strict).
% minScoreRatioG : double
%     a gene is only assigned to KOs for which the score is
%     >=log(score)/log(best score) for that gene. This is to prevent that a
%     gene which clearly belongs to one KO is assigned also to KOs with much
%     lower scores (default 0.8, lower is less strict).
% maxPhylDist : double
%     -1 to only use sequences from the same domain (Prokaryota, Eukaryota);
%     any other (positive) value to only use sequences for organisms where
%     the phylogenetic distance is at the most this large (as calculated in
%     getPhylDist). Only used when reconstructing from KEGG annotations (no
%     fastaFile) (default Inf, which means that all sequences will be used).
% globalModel : struct
%     structure containing both model and KOModel structures as generated
%     by getModelFromKEGG. These will otherwise be loaded by via
%     getModelFromKEGG. Providing globalKEGGmodel can speed up model
%     generation if getKEGGModelForOrganism is run multiple times for
%     different strains. Example:
%     [globalModel.model,globalModel.KOModel] = getModelFromKEGG (default
%     empty, global model is loaded by getModelFromKEGG).
%
% Returns
% -------
% model : struct
%     the reconstructed model.
%
% Notes
% -----
% The reconstruction works in one of two modes:
%
%   1. From KEGG annotations (no fastaFile supplied). The reactions
%      associated with organismID in the local KEGG database are kept;
%      maxPhylDist controls which organisms' annotations are considered.
%   2. From protein homology (fastaFile supplied). The query proteome is
%      searched, in a single hmmsearch, against a prebuilt KEGG-version- and
%      domain-specific concatenated KO HMM library (e.g. kegg118_eukaryotes),
%      downloaded from the corresponding RAVEN release if not already present
%      in dataDir. Hits are filtered by cutOff and the minScoreRatioKO /
%      minScoreRatioG ratios into a KO-gene matrix, from which the model is
%      built (spontaneous reactions without genes are kept if keepSpontaneous
%      is true).
%
% The global KEGG model (reaction-KO relationships) is either supplied via
% globalModel, built from a local KEGG FTP dump in dataDir\keggdb (which
% requires a KEGG FTP license), or loaded from the internal RAVEN version.
%
% Examples
% --------
%     model=getKEGGModelForOrganism(organismID,fastaFile,dataDir,...
%        outDir,keepSpontaneous,keepUndefinedStoich,keepIncomplete,...
%        keepGeneral,cutOff,minScoreRatioKO,minScoreRatioG,maxPhylDist);

p=parseRAVENargs(varargin, {'fastaFile',[]; 'dataDir',[]; 'outDir',[]; ...
    'keepSpontaneous',true; 'keepUndefinedStoich',true; ...
    'keepIncomplete',true; 'keepGeneral',false; 'cutOff',10^-50; ...
    'minScoreRatioKO',0.3; 'minScoreRatioG',0.8; 'maxPhylDist',inf; ...
    'globalModel',[]});
fastaFile=p.fastaFile;
if isempty(fastaFile)
    fastaFile=[];
else
    fastaFile=char(fastaFile);
end
dataDir=p.dataDir;
if isempty(dataDir)
    dataDir=[];
else
    dataDir=char(dataDir);
end
outDir=p.outDir;
if isempty(outDir)
    outDir=tempdir;
    %Delete all *.out files if any exist
    delete(fullfile(outDir,'*.out'));
else
    outDir=char(outDir);
end
keepSpontaneous=p.keepSpontaneous;
if isempty(keepSpontaneous)
    keepSpontaneous=true;
end
keepUndefinedStoich=p.keepUndefinedStoich;
if isempty(keepUndefinedStoich)
    keepUndefinedStoich=true;
end
keepIncomplete=p.keepIncomplete;
if isempty(keepIncomplete)
    keepIncomplete=true;
end
keepGeneral=p.keepGeneral;
if isempty(keepGeneral)
    keepGeneral=false;
end
cutOff=p.cutOff;
if isempty(cutOff)
    cutOff=10^-50;
end
minScoreRatioKO=p.minScoreRatioKO;
if isempty(minScoreRatioKO)
    minScoreRatioKO=0.3;
end
minScoreRatioG=p.minScoreRatioG;
if isempty(minScoreRatioG)
    minScoreRatioG=0.8;
end
maxPhylDist=p.maxPhylDist;
if isempty(maxPhylDist)
    maxPhylDist=inf;
    %Include all sequences for each reaction
end
globalModel=p.globalModel;

if isempty(fastaFile)
    fprintf(['\n*** The model reconstruction from KEGG based on the annotation available for KEGG Species <strong>' organismID '</strong> ***\n\n']);
else
    fprintf('\n*** The model reconstruction from KEGG based on the protein homology search against KEGG Orthology specific HMMs ***\n\n');
    %Check if query fasta exists
    fastaFile=checkFileExistence(fastaFile,2); %Copy file to temp dir
end

%Run the external binaries multi-threaded to use all logical cores assigned
%to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Get the directory for RAVEN Toolbox.
ravenPath=findRAVENroot();

libraryFile='';

%Checking if dataDir is consistent. It must point to a RAVEN-provided
%pre-trained HMM set, compatible with the current RAVEN version. The
%reconstruction uses the concatenated KO HMM library (a single
%gzip-compressed flatfile, queried in one hmmsearch); if it is not already
%present it is downloaded from the corresponding raven-data release
%(https://github.com/SysBioChalmers/raven-data).
if ~isempty(dataDir)
    hmmOptions={'kegg118_eukaryotes','kegg118_prokaryotes'};
    if ~endsWith(dataDir,hmmOptions) %Check if dataDir ends with any of the hmmOptions.
        %If not, then check whether the required keggdb folder exists anyway.
        if ~isfile(fullfile(dataDir,'keggdb','genes.pep'))
            error(['Pre-trained HMMs set is not recognised. If you want to download RAVEN provided sets, it should match any of the following: ' strjoin(hmmOptions,' or ')])
        end
    else
        %dataDir points to a RAVEN-provided set. Use the concatenated KO HMM
        %library (one gzip-compressed flatfile, queried in a single
        %hmmsearch), downloading and extracting it if necessary. The dataDir
        %name matches the published HMM library asset, so it doubles as the
        %download filename.
        hmmName=hmmOptions{endsWith(dataDir,hmmOptions)};
        libraryFile=[dataDir '.hmm'];
        if isfile(libraryFile)
            fprintf(['NOTE: Found <strong>' libraryFile '</strong> HMM library, it will therefore be used during reconstruction\n']);
        elseif isfile([libraryFile '.gz'])
            fprintf('Extracting the HMM library file... ');
            gunzip([libraryFile '.gz']);
            fprintf('COMPLETE\n');
            useConcatLib=false;
        else
            fprintf('Downloading the HMM library file... ');
            try
                websave([libraryFile '.gz'],['https://github.com/SysBioChalmers/raven-data/releases/download/kegg118/' hmmName '.hmm.gz']);
            catch ME
                if strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError')
                    error('Failed to download the HMM library file, the server returned a 404 error, try again later. If the problem persists please report it on the RAVEN GitHub Issues page: https://github.com/SysBioChalmers/RAVEN/issues')
                end
            end
            fprintf('COMPLETE\n');
            fprintf('Extracting the HMM library file... ');
            gunzip([libraryFile '.gz']);
            fprintf('COMPLETE\n');
            useConcatLib=true;
        end
        %Check that the HMM library is available
        if ~isfile(libraryFile)
            error(['The HMM library seems improperly extracted and not found at ',strrep(libraryFile,'\','/'),'. Please remove the corresponding .gz file and rerun getKEGGModelForOrganism']);
        end
    end
end

%Check if the fasta-file contains '/' or'\'. If not then it's probably just
%a file name. Expand to full path.
if any(fastaFile)
    if ~any(strfind(fastaFile,'\')) && ~any(strfind(fastaFile,'/'))
        fastaFile=which(fastaFile);
    end
    %Create the required sub-folders in dataDir if they dont exist
    if ~isfolder(fullfile(dataDir,'keggdb'))
        mkdir(dataDir,'keggdb');
    end
    if ~isfolder(outDir)
        mkdir(outDir);
    end
end

%First generate the full global KEGG model. Can be provided as input.
%Otherwise, getModelFromKEGG is run. The dataDir must not be supplied as
%there is also an internal RAVEN version available
if ~isempty(globalModel)
    model=globalModel.model;
    KOModel=globalModel.KOModel;
elseif any(dataDir)
    [model, KOModel]=getModelFromKEGG(fullfile(dataDir,'keggdb'),keepSpontaneous,keepUndefinedStoich,keepIncomplete,keepGeneral);
else
    [model, KOModel]=getModelFromKEGG([],keepSpontaneous,keepUndefinedStoich,keepIncomplete,keepGeneral);
end
model.id=organismID;
model.c=zeros(numel(model.rxns),1);

%If no FASTA file is supplied, then just remove all genes which are not for
%the given organism ID
if isempty(fastaFile)
    %Check if organismID can be found in KEGG species list or is
    %set to "eukaryotes" or "prokaryotes"
    phylDistsFull=getPhylDist(fullfile(dataDir,'keggdb'),true);
    if ~ismember(organismID,[phylDistsFull.ids 'eukaryotes' 'prokaryotes'])
        error('Provided organismID is incorrect. Only species abbreviations from KEGG Species List or "eukaryotes"/"prokaryotes" are allowed.');
    end

    fprintf(['Pruning the model from <strong>non-' organismID '</strong> genes... ']);
    if ismember(organismID,{'eukaryotes','prokaryotes'})
        phylDists=getPhylDist(fullfile(dataDir,'keggdb'),maxPhylDist==-1);
        if strcmp(organismID,'eukaryotes')
            proxyid='hsa';
            %Use H. sapiens here
        else
            proxyid='eco';
            %Use E. coli here
        end
        [~, phylDistId]=ismember(proxyid,phylDists.ids);
        idsToKeep=phylDists.ids(~isinf(phylDists.distMat(phylDistId,:)));
        taxIDs=cellfun(@(x) x{1},cellfun(@(x) strsplit(x,':'),model.genes,'UniformOutput',false),'UniformOutput',false);
        I=ismember(upper(taxIDs),upper(idsToKeep));
    else
        %KEGG organism IDs may have three or four letters
        organismID=strcat(organismID,':');
        %Add colon for accurate matching
        if length(organismID)==4
            I=cellfun(@(x) strcmpi(x(1:4),organismID),model.genes);
        elseif length(organismID)==5
            I=cellfun(@(x) strcmpi(x(1:5),organismID),model.genes);
        end
    end
    %Remove those genes
    model.genes=model.genes(I);
    model.rxnGeneMat=model.rxnGeneMat(:,I);
    fprintf('COMPLETE\n');
end

%First remove all reactions without genes
if keepSpontaneous==true
    fprintf('Removing non-spontaneous reactions without GPR rules... ');
    load(fullfile(ravenPath,'reconstruction','kegg','keggRxns.mat'),'isSpontaneous');
    I=~any(model.rxnGeneMat,2)&~ismember(model.rxns,isSpontaneous);
    spontRxnsWithGenes=model.rxns(any(model.rxnGeneMat,2)&~ismember(model.rxns,isSpontaneous));
else
    fprintf('Removing reactions without GPR rules... ');
    I=~any(model.rxnGeneMat,2);
end
model=removeReactions(model,I,true);
fprintf('COMPLETE\n');

%Clean gene names
fprintf('Fixing gene names in the model... ');
%Get rid of the prefix organism id
model.genes=regexprep(model.genes,'^\w+?:','');
fprintf('COMPLETE\n');

%If no FASTA file is supplied, then we are done here
if isempty(fastaFile)
    %Create grRules
    fprintf('Constructing GPR associations and annotations for the model... ');
    model.grRules=cell(numel(model.rxns),1);
    model.grRules(:)={''};
    %Add the gene associations as 'or'
    for i=1:numel(model.rxns)
        %Find the involved genes
        I=find(model.rxnGeneMat(i,:));
        if any(I)
            model.grRules{i}=['(' model.genes{I(1)}];
            for j=2:numel(I)
                model.grRules{i}=[model.grRules{i} ' or ' model.genes{I(j)}];
            end
            model.grRules{i}=[model.grRules{i} ')'];
        end
    end
    %Fix grRules and reconstruct rxnGeneMat
    [grRules,rxnGeneMat] = standardizeGrRules(model); %Give detailed output
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;
    %Add geneMiriams, assuming that it follows the syntax
    %kegg.genes/organismID:geneName
    model.geneMiriams='';
    for i=1:numel(model.genes)
        model.geneMiriams{i,1}.name{1,1}='kegg.genes';
        model.geneMiriams{i,1}.value{1,1}=strcat(lower(organismID),model.genes{i,1});
    end
    %Add the description to the reactions
    for i=1:numel(model.rxns)
        if ~isempty(model.rxnNotes{i})
            model.rxnNotes(i)=strcat('Included by getKEGGModelForOrganism (without HMMs).',model.rxnNotes(i));
            model.rxnNotes(i)=strrep(model.rxnNotes(i),'.','. ');
        else
            model.rxnNotes(i)={'Included by getKEGGModelForOrganism (without HMMs)'};
        end
    end
    fprintf('COMPLETE\n\n');
    fprintf('*** Model reconstruction complete ***\n');
    return;
end

%Query the whole proteome against the concatenated KO HMM library in a
%single hmmsearch. With the profile library as the query and the
%proteome as the target sequence database, the reported per-hit
%E-values match RAVEN's historical per-KO hmmsearch (same search
%direction, same effective database size), so the downstream scoring is
%unchanged - thousands of hmmsearch calls simply collapse into one, and
%no per-organism phylogenetic-distance subsampling is needed because the
%library sequence set is already fixed.
if ismac
    binEnd='.mac';
else
    binEnd='';
end
tblFile=[tempname '.tblout'];
fprintf('Querying the user-specified FASTA file against the KEGG Orthology specific HMMs... ');
if ismac || isunix
    [status, output]=system(['"' fullfile(ravenPath,'software','hmmer',['hmmsearch' binEnd]) '" --cpu "' num2str(cores) '" --tblout "' tblFile '" "' libraryFile '" "' fastaFile '"']);
else
    wslPath.hmmsearch   = getWSLpath(fullfile(ravenPath,'software','hmmer','hmmsearch'));
    wslPath.libraryFile = getWSLpath(libraryFile);
    wslPath.fastaFile   = getWSLpath(fastaFile);
    wslPath.tblFile     = getWSLpath(tblFile);
    [status, output]=system(['wsl "' wslPath.hmmsearch '" --cpu "' num2str(cores) '" --tblout "' wslPath.tblFile '" "' wslPath.libraryFile '" "' wslPath.fastaFile '"']);
end
if status~=0
    EM=['Error when querying the concatenated HMM library:\n' output];
    dispEM(EM);
end
fprintf('COMPLETE\n');

fprintf('Parsing the HMM search results... ');
koGeneMat=zeros(numel(KOModel.rxns),3000); %Make room for 3000 genes
genes=cell(3000,1);
%Map each KO id to its row in koGeneMat
koTable=java.util.Hashtable;
for i=1:numel(KOModel.rxns)
    koTable.put(KOModel.rxns{i},i);
end
%Store the column index for each gene in a hash list
hTable=java.util.Hashtable;
geneCounter=0;
fid=fopen(tblFile,'r');
while 1
    tline=fgetl(fid);
    if ~ischar(tline)
        break;
    end
    %Skip the comment/header lines that hmmsearch --tblout writes
    if isempty(tline) || tline(1)=='#'
        continue;
    end
    elements=regexp(strtrim(tline),'\s+','split');
    if numel(elements)<5
        continue;
    end
    %With the profile library as the query, the target name (col 1) is
    %the proteome gene, the query name (col 3) is the KO, and the
    %full-sequence E-value is in col 5
    gene=elements{1};
    ko=elements{3};
    score=str2double(elements{5});
    koIdx=koTable.get(ko);
    if isempty(koIdx) || isnan(score) || score>cutOff
        continue;
    end
    %If the score is exactly 0, change it to a very small value to avoid
    %NaN
    if score==0
        score=10^-250;
    end
    I=hTable.get(gene);
    if any(I)
        %Keep the best (smallest) E-value for this gene-KO pair
        if koGeneMat(koIdx,I)==0 || score<koGeneMat(koIdx,I)
            koGeneMat(koIdx,I)=score;
        end
    else
        geneCounter=geneCounter+1;
        hTable.put(gene,geneCounter);
        genes{geneCounter}=gene;
        koGeneMat(koIdx,geneCounter)=score;
    end
end
fclose(fid);
delete(tblFile);
fprintf('COMPLETE\n');
end

fprintf('Removing gene, KEGG Orthology associations below minScoreRatioKO, minScoreRatioG... ');
koGeneMat=koGeneMat(:,1:geneCounter);

%Remove the genes for each KO that are below minScoreRatioKO.
for i=1:size(koGeneMat,1)
    J=find(koGeneMat(i,:));
    if any(J)
        koGeneMat(i,J(log(koGeneMat(i,J))/log(min(koGeneMat(i,J)))<minScoreRatioKO))=0;
    end
end

%Remove the KOs for each gene that are below minScoreRatioG
for i=1:size(koGeneMat,2)
    J=find(koGeneMat(:,i));
    if any(J)
        koGeneMat(J(log(koGeneMat(J,i))/log(min(koGeneMat(J,i)))<minScoreRatioG),i)=0;
    end
end
fprintf('COMPLETE\n');

fprintf('Adding gene annotations to the model... ');
%Create the new model
model.genes=genes(1:geneCounter);
model.grRules=cell(numel(model.rxns),1);
model.grRules(:)={''};
model.rxnGeneMat=sparse(numel(model.rxns),numel(model.genes));

%Loop through the reactions and add the corresponding genes
for i=1:numel(model.rxns)
    if isstruct(model.rxnMiriams{i})
        %Get all KOs
        I=find(strcmpi(model.rxnMiriams{i}.name,'kegg.orthology'));
        KOs=model.rxnMiriams{i}.value(I);
        %Find the KOs and the corresponding genes
        J=ismember(KOModel.rxns,KOs);
        [~, K]=find(koGeneMat(J,:));

        if any(K)
            model.rxnGeneMat(i,K)=1;
            %Also delete KOs for which no genes were found. If no genes at
            %all were matched to the reaction it will be deleted later
            L=sum(koGeneMat(J,:),2)==0;
            model.rxnMiriams{i}.value(I(L))=[];
            model.rxnMiriams{i}.name(I(L))=[];
        end
    end
end
fprintf('COMPLETE\n');

%Find and delete all reactions without genes. This also removes genes that
%are not used (which could happen because minScoreRatioG and
%minScoreRatioKO). If keepSpontaneous==true, the spontaneous reactions
%without genes are kept in the model. Spontaneous reactions with original
%gene associations are treated in the same way, like the rest of the
%reactions - if gene associations were removed during HMM search, such
%reactions are deleted from the model
if keepSpontaneous==true
    %Not the most comprise way to delete reactions without genes, but this
    %makes the code easier to understand. Firstly the non-spontaneous
    %reactions without genes are removed. After that, the second deletion
    %step removes spontaneous reactions, which had gene associations before
    %HMM search, but no longer have after it
    fprintf('Removing non-spontaneous reactions which after HMM search no longer have GPR rules... ');
    I=~any(model.rxnGeneMat,2)&~ismember(model.rxns,isSpontaneous);
    model=removeReactions(model,I,true,true);
    I=~any(model.rxnGeneMat,2)&ismember(model.rxns,spontRxnsWithGenes);
    model=removeReactions(model,I,true,true);
else
    %Just simply check for any new reactions without genes and remove
    %it
    fprintf('Removing reactions which after HMM search no longer have GPR rules... ');
    I=~any(model.rxnGeneMat,2);
    model=removeReactions(model,I,true,true);
end
fprintf('COMPLETE\n');

fprintf('Constructing GPR rules and finalizing the model... ');
%Add the gene associations as 'or'
for i=1:numel(model.rxns)
    %Find the involved genes
    I=find(model.rxnGeneMat(i,:));
    if any(I)
        model.grRules{i}=['(' model.genes{I(1)}];
        for j=2:numel(I)
            model.grRules{i}=[model.grRules{i} ' or ' model.genes{I(j)}];
        end
        model.grRules{i}=[model.grRules{i} ')'];
    end
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,false); %Give detailed output
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;

%Add the description to the reactions
for i=1:numel(model.rxns)
    if ~isempty(model.rxnNotes{i})
        model.rxnNotes(i)=strcat('Included by getKEGGModelForOrganism (using HMMs).',model.rxnNotes(i));
        model.rxnNotes(i)=strrep(model.rxnNotes(i),'.','. ');
    else
        model.rxnNotes(i)={'Included by getKEGGModelForOrganism (using HMMs)'};
    end
end
%Remove the temp fasta file
delete(fastaFile)
fprintf('COMPLETE\n\n*** Model reconstruction complete ***\n');
end

function files=listFiles(directory)
%Supporter function to list the files in a directory and return them as a
%cell array
temp=dir(directory);
files=cell(numel(temp),1);
for i=1:numel(temp)
    files{i}=temp(i,1).name;
end
files=strrep(files,'.fa','');
files=strrep(files,'.hmm','');
files=strrep(files,'.out','');
files=strrep(files,'.faw','');
end
