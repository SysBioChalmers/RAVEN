function model=getKEGGModelForOrganism(organismID,fastaFile,dataDir,outDir,...
    keepUndefinedStoich,keepIncomplete,keepGeneral,cutOff,minScoreRatioG,...
    minScoreRatioKO,maxPhylDist,nSequences,seqIdentity)
% getKEGGModelForOrganism
%   Reconstructs a genome-scale metabolic model based on protein homology to the
%   orthologies in KEGG
%
%   organismID          three or four letter abbreviation of the organism
%                       (as used in KEGG). If not available, use a closely
%                       related species. This is used for determing the
%                       phylogenetic distance. Use 'eukaryotes' or
%                       'prokaryotes' to get a model for the whole domain.
%                       Only applicable if fastaFile is empty, i.e. no
%                       homology search should be performed
%   fastaFile           a FASTA file that contains the protein sequences of
%                       the organism for which to reconstruct a model (opt,
%                       if no FASTA file is supplied then a model is
%                       reconstructed based only on the organism
%                       abbreviation. This option ignores all settings
%                       except for keepUndefinedStoich, keepIncomplete and
%                       keepGeneral)
%   dataDir             directory for which to retrieve the input data.
%                       Should contain a combination of these sub-folders:
%                       -dataDir\keggdb
%                           The KEGG database files used in 1a (see below)
%                       -dataDir\fasta
%                           The multi-FASTA files generated in 1b (see
%                           below)
%                       -dataDir\aligned
%                           The aligned FASTA files as generated in 2a (see
%                           below)
%                       -dataDir\hmms
%                           The hidden Markov models as generated in 2b or
%                           downloaded from BioMet Toolbox (see below)
%                       The final directory in dataDir should be styled as
%                       proXXX_keggYY or eukXXX_keggYY, indicating whether
%                       the HMMs were trained on pro- or eukaryotic
%                       sequences, using a sequence similarity threshold of
%                       XXX %, fitting the KEGG version YY. E.g.
%                       euk100_kegg82. (opt, see note about fastaFile. Note
%                       that in order to rebuild the KEGG model from a
%                       database dump, as opposed to using the version
%                       supplied with RAVEN, you would still need to supply
%                       this)
%   outDir              directory to save the results from the quering of
%                       the hidden Markov models. The output is specific
%                       for the input sequences and the settings used. It
%                       is stored in this manner so that the function can
%                       continue if interrupted or if it should run in
%                       parallel. Be careful not to leave output files from
%                       different organisms or runs with different settings
%                       in the same folder. They will not be overwritten
%                       (opt, default is a temporary dir where all *.out
%                       files are deleted before and after doing the
%                       reconstruction)
%   keepUndefinedStoich	include reactions in the form n A <=> n+1 A. These
%                       will be dealt with as two separate metabolites
%                       (opt, default true)
%   keepIncomplete      include reactions which have been labelled as
%                       "incomplete", "erroneous" or "unclear" (opt,
%                       default true)
%   keepGeneral         include reactions which have been labelled as
%                       "general reaction". These are reactions on the form
%                       "an aldehyde <=> an alcohol", and are therefore
%                       unsuited for modelling purposes. Note that not all
%                       reactions have this type of annotation, and the
%                       script will therefore not be able to remove all
%                       such reactions (opt, default false)
%   cutOff              significance score from HMMer needed to assign
%                       genes to a KO (opt, default 10^-50)
%   minScoreRatioG      a gene is only assigned to KOs for which the score
%                       is >=log(score)/log(best score) for that gene. This
%                       is to prevent thata gene which clearly belongs to
%                       one KO is assigned also to KOs with much lower
%                       scores (opt, default 0.8 (lower is less strict))
%   minScoreRatioKO     ignore genes in a KO if their score is
%                       <log(score)/log(best score in KO). This is to
%                       "prune" KOs which havemany genes and where some are
%                       clearly a better fit (opt, default 0.3 (lower is
%                       less strict))
%   maxPhylDist         -1: only use sequences from the same domain
%                       (Prokaryota, Eukaryota)
%                       other (positive) value: only use sequences for
%                       organisms where the phylogenetic distance is at the
%                       most this large (as calculated in getPhylDist)
%                       (opt, default Inf, which means that all sequences
%                       will be used)
%   nSequences          for each KO, use up to this many sequences from the
%                       most closely related species. This is mainly to
%                       speed up the alignment process for KOs with very
%                       many genes. This subsampling is performed before
%                       running CD-HIT (opt, default inf)
%   seqIdentity         sequence identity threshold in CD-HIT, referred as
%                       "global sequence identity" in CD-HIT User's Guide.
%                       The only possible options are 1 (100 %), 0.9 (90 %)
%                       and 0.5 (50 %). If other values are provided,
%                       CD-HIT is skipped (opt, default -1, i.e. CD-HIT is
%                       skipped)
%
%   model               the reconstructed model
%
%   PLEASE READ THIS: The input to this function can be confusing, because
%   it is intended to be run in parallel on a cluster or in multiple
%   sessions. It therefore saves a lot of intermediate results to storage.
%   This also serves the purpose of not having to do redundant
%   calculations. This, however, comes with the disadvantage of somewhat
%   trickier handling. This is what this function does:
%
%   1a. Loads files from a local KEGG FTP dump and constructs a general
%       RAVEN model representing the metabolic network. KEGG FTP access
%       requires a <a href="matlab:
%       web('http://www.bioinformatics.jp/en/keggftp.html')">license</a>.
%   1b. Generates FASTA files from the KEGG FTP dump (see 1a). One
%       multi-FASTA file for each KO in KEGG is generated. Make sure you
%       actually need these first, as the parsing is time-consuming (see
%       below).
%
%   These steps only have to be re-done every time KEGG updates their
%   database (or rather when the updates are large enough to warrant
%   rerunning this part). Many user would probably never use this feature.
%
%   2a. Does alignment of the multi-FASTA files for future use. This uses
%       the settings "useEvDist" and "nSequences" to control which
%       sequences should be used for constructing Hidden Markov models
%       (HMMs), and later for matching your sequences to.
%   2b. Trains hidden Markov models using HMMer for each of the aligned
%       FASTA files.
%
%   The most common alternatives here would be to use sequences from only
%   eukaryotes, only prokaryotes or all sequences in KEGG. As explained in
%   the README.md file, various sets of pre-trained hidden Markov models are
%   available at <a href="matlab:
%   web('http://biomet-toolbox.chalmers.se/index.php?page=downtools-raven')">BioMet
%   Toolbox</a>. This is normally the most convenient way, but if you would
%   like to use, for example, only fungal sequences for training the HMMs
%   then you need to run this part.
%
%   3a. Queries the HMMs with sequences for the organism you are making a
%       model for. This step uses both the output from step 1a and from 2b.
%   3b. Constructs a model based on the fit to the HMMs and the chosen
%       parameters.
%
%   These steps are specific to the organism for which you are
%   reconstructing the model.
%
%   In principle the function looks at which output that is already available
%   and runs only the parts that are required for step 3. This means
%   that (see the definition of the parameters for details):
%   -1a is only performed if there are no KEGG model files in the
%   RAVEN\external\kegg directory
%   -1b is only performed if not all required HMMs OR aligned FASTA files
%   OR multi-FASTA files exist in the defined dataDir. This means that this
%   step is skipped if the HMMs are downloaded from BioMet Toolbox instead
%   (see above). If not all files exist it will try to find
%   the KEGG database files in dataDir.
%   -2a is only performed if not all required HMMs OR aligned FASTA files
%   files exist in the defined dataDir. This means that this step is skipped
%   if the HMMs are downloaded from BioMet Toolbox instead (see above).
%   -2b is only performed if not all required HMMs exist in the defined
%   dataDir. This means that this step is skipped if the FASTA files or
%   HMMs are downloaded from BioMet Toolbox instead (see above).
%   -3a is performed for the required HMMs for which no corresponding .out
%   file exists in outDir. This is just a way to enable the function to be
%   run in parallel or to resume if interrupted.
%   -3b is always performed.
%
%   Usage: model=getKEGGModelForOrganism(organismID,fastaFile,dataDir,outDir,...
%    keepUndefinedStoich,keepIncomplete,keepGeneral,cutOff,minScoreRatioG,...
%    minScoreRatioKO,maxPhylDist,nSequences,seqIdentity)
%
%   Simonas Marcisauskas, 2018-04-03
%

if nargin<2
    fastaFile=[];
end
if nargin<3
    dataDir=[];
end
if nargin<4
    outDir=[];
end
if isempty(outDir)
    outDir=tempdir;
    %Delete all *.out files if any exist
    delete(fullfile(outDir,'*.out'));
end
if nargin<5
    keepUndefinedStoich=true;
end
if nargin<6
    keepIncomplete=true;
end
if nargin<7
    keepGeneral=false;
end
if nargin<8
    cutOff=10^-50;
end
if nargin<9
    minScoreRatioG=0.8;
end
if nargin<10
    minScoreRatioKO=0.3;
end
if nargin<11
    maxPhylDist=inf;
    %Include all sequences for each reaction
end
if nargin<12
    nSequences=inf;
    %Include all sequences for each reaction
end
if nargin<13
    seqIdentity=-1;
    %CD-HIT is not used in the pipeline
end

%Run the external binaries multi-threaded to use all logical cores assigned
%to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Checking if dataDir is consistent. It must point to pre-trained HMMs set,
%compatible with the the current RAVEN version.
%The user may have the required zip file already in working directory or
%have it extracted. If the zip file and directory is not here, it is
%downloaded from BioMet ToolBox 2.0 server.
if ~isempty(dataDir)
    hmmOptions={'euk100_kegg82'; ...
        'euk90_kegg82'; ...
        'euk50_kegg82'; ...
        'prok100_kegg82'; ...
        'prok90_kegg82'; ...
        'prok50_kegg82'};
    if all(cellfun(@isempty,regexp(dataDir,strcat(hmmOptions,'$')))) % Check if dataDir ends with any of the hmmOptions
        EM='Pre-trained HMMs set is not recognised. It should match any of the following sets (which are available to download):';
        disp(EM);
        disp(hmmOptions);
        return;
    end
    if exist(dataDir,'dir')
        fprintf('Provided dataDir is in correct format for this RAVEN version in order to use pre-trained HMMs...\n')
    elseif ~exist(dataDir,'dir') && exist([dataDir,'.zip'],'file')
        fprintf('Extracting HMMs archive file...\n');
        unzip([dataDir,'.zip']);
    else
        fprintf('Downloading HMMs archive file...\n');
        websave([dataDir,'.zip'],['http://biomet-toolbox.chalmers.se/tools/downloadable/files/',dataDir,'.zip']);
        fprintf('Extracting HMMs archive file...\n');
        unzip([dataDir,'.zip']);
    end
end

%Check if the fasta-file contains '/' or'\'. If not then it's probably just
%a file name. It is then merged with the current folder
if any(fastaFile)
    if ~any(strfind(fastaFile,'\')) && ~any(strfind(fastaFile,'/'))
        fastaFile=fullfile(pwd,fastaFile);
    end
    %Create the required sub-folders in dataDir if they dont exist
    if ~exist(fullfile(dataDir,'keggdb'),'dir')
        mkdir(dataDir,'keggdb');
    end
    if ~exist(fullfile(dataDir,'fasta'),'dir')
        mkdir(dataDir,'fasta');
    end
    if ~exist(fullfile(dataDir,'aligned'),'dir')
        mkdir(dataDir,'aligned');
    end
    if ~exist(fullfile(dataDir,'hmms'),'dir')
        mkdir(dataDir,'hmms');
    end
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
end

%First generate the full KEGG model. The dataDir mustn't be supplied as
%there is also an internal RAVEN version available
if any(dataDir)
    [model, KOModel]=getModelFromKEGG(fullfile(dataDir,'keggdb'),keepUndefinedStoich,keepIncomplete,keepGeneral);
else
    [model, KOModel]=getModelFromKEGG([],keepUndefinedStoich,keepIncomplete,keepGeneral);
end
fprintf('Completed generation of KEGG model\n');
model.id=organismID;
model.c=zeros(numel(model.rxns),1);

%If no FASTA file is supplied, then just remove all genes which are not for
%the given organism ID
if isempty(fastaFile)
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
end

%First remove all reactions without genes
hasGenes=any(model.rxnGeneMat,2);
model=removeReactions(model,~hasGenes,true);

%Clean gene names
for i=1:numel(model.genes)
    %First get rid of the prefix organism id
    model.genes{i}=model.genes{i}(strfind(model.genes{i},':')+1:end);
    %Find and remove the description in parentheses if any
    s=strfind(model.genes{i},'(');
    if any(s)
        model.genes{i}=model.genes{i}(1:s-1);
    end
end

%If no FASTA file is supplied, then we're done here
if isempty(fastaFile)
    %Create grRules
    model.grRules=cell(numel(model.rxns),1);
    model.grRules(:)={''};
    %Add the gene associations as 'or'
    for i=1:numel(model.rxns)
        %Find the involved genes
        I=find(model.rxnGeneMat(i,:));
        model.grRules{i}=['(' model.genes{I(1)}];
        for j=2:numel(I)
            model.grRules{i}=[model.grRules{i} ' or ' model.genes{I(j)}];
        end
        model.grRules{i}=[model.grRules{i} ')'];
    end
    %Add geneMiriams, assuming that it follows the syntax
    %kegg.genes/organismID:geneName
    model.geneMiriams='';
    for i=1:numel(model.genes)
        model.geneMiriams{i,1}.name{1,1}='kegg.genes';
        model.geneMiriams{i,1}.value{1,1}=strcat(lower(organismID),model.genes{i,1});
    end
    return;
end

%Trim the genes so that they only contain information that can be matched
%to the KEGG file of protein sequences (remove all information after first
%parenthesis)
%NOTE: For some reason the organism abbreviation should be in lower case in
%this database. Fix this here
for i=1:numel(KOModel.genes)
    parIndex=strfind(KOModel.genes{i},'(');
    if any(parIndex)
        KOModel.genes{i}=KOModel.genes{i}(1:parIndex-1);
    end
    colIndex=strfind(KOModel.genes{i},':');
    KOModel.genes{i}=[lower(KOModel.genes{i}(1:colIndex-1)) KOModel.genes{i}(colIndex:end)];
end

%Create a phylogenetic distance structure
phylDistStruct=getPhylDist(fullfile(dataDir,'keggdb'),maxPhylDist==-1);
[~, phylDistId]=ismember(model.id,phylDistStruct.ids);
fprintf('Completed creation of phylogenetic distance matrix\n');

%Calculate the real maximal distance now. An abitary large number of 1000
%is used for the "all in kingdom" or "all sequences" options. This is a bit
%inconvenient way to do it, but it's to make it fit with some older code
if isinf(maxPhylDist) || maxPhylDist==-1
    maxPhylDist=1000;
end

%Get the KO ids for which files have been generated. Maybe not the neatest
%way..
fastaFiles=listFiles(fullfile(dataDir,'fasta','*.fa'));
alignedFiles=listFiles(fullfile(dataDir,'aligned','*.fa'));
alignedWorking=listFiles(fullfile(dataDir,'aligned','*.faw'));
hmmFiles=listFiles(fullfile(dataDir,'hmms','*.hmm'));
outFiles=listFiles(fullfile(outDir,'*.out'));

%Check if multi-FASTA files should be generated. This should only be
%performed if there are IDs in the KOModel structure that haven't been
%parsed yet
missingFASTA=setdiff(KOModel.rxns,[fastaFiles;alignedFiles;hmmFiles;outFiles]);

if ~isempty(missingFASTA)
    if ~exist(fullfile(dataDir,'keggdb','genes.pep'),'file')
        EM=fprintf(['The file ''genes.pep'' cannot be located at ' strrep(dataDir,'\','/') '/ and should be downloaded from the KEGG FTP.\n']);
        dispEM(EM);
    end
    %Only construct models for KOs which don't have files already
    fastaModel=removeReactions(KOModel,setdiff(KOModel.rxns,missingFASTA),true,true);
    %Permute the order of the KOs in the model so that constructMultiFasta
    %can be run on several processors at once
    fastaModel=permuteModel(fastaModel,randperm(RandStream.create('mrg32k3a','Seed',cputime()),numel(fastaModel.rxns)),'rxns');
    constructMultiFasta(fastaModel,fullfile(dataDir,'keggdb','genes.pep'),fullfile(dataDir,'fasta'));
end
fprintf('Completed generation of multi-FASTA files\n');

%Get the directory for RAVEN Toolbox. This is to get the path to the third
%party software used
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='';
else
    EM='Unknown OS, exiting.';
    disp(EM);
    return
end

%Check if alignment of FASTA files should be performed
missingAligned=setdiff(KOModel.rxns,[alignedFiles;hmmFiles;alignedWorking;outFiles]);
if ~isempty(missingAligned)
    missingAligned=missingAligned(randperm(RandStream.create('mrg32k3a','Seed',cputime()),numel(missingAligned)));
    %Align all sequences using MAFFT
    for i=1:numel(missingAligned)
        %This is checked here because it could be that it is created by a
        %parallel process. The faw-files are saved as temporary files to
        %kept track of which files are being worked on
        if ~exist(fullfile(dataDir,'aligned',[missingAligned{i} '.faw']),'file') && ~exist(fullfile(dataDir,'aligned',[missingAligned{i} '.fa']),'file')
            %Check that the multi-FASTA file exists. It should do so since
            %we are saving empty files as well. Print a warning and
            %continue if not
            if ~exist(fullfile(dataDir,'fasta',[missingAligned{i} '.fa']),'file')
                EM=['WARNING: The multi-FASTA file for ' missingAligned{i} ' does not exist'];
                dispEM(EM,false);
                continue;
            end
            
            %If the multi-FASTA file is empty then save an empty aligned
            %file and continue
            s=dir(fullfile(dataDir,'fasta',[missingAligned{i} '.fa']));
            if s.bytes<=0
                fid=fopen(fullfile(dataDir,'aligned',[missingAligned{i} '.fa']),'w');
                fclose(fid);
                continue;
            end
            
            %Create an empty file to prevent other threads to start to work on the same alignment
            fid=fopen(fullfile(dataDir,'aligned',[missingAligned{i} '.faw']),'w');
            fclose(fid);
            
            %First load the FASTA file, then select up to nSequences sequences
            %of the most closely related species, apply any constraints from
            %maxPhylDist, and save it as a temporary file, and create the
            %model from that
            
            fastaStruct=fastaread(fullfile(dataDir,'fasta',[missingAligned{i} '.fa']));
            phylDist=inf(numel(fastaStruct),1);
            for j=1:numel(fastaStruct)
                %Get the organism abbreviation
                index=strfind(fastaStruct(j).Header,':');
                if any(index)
                    abbrev=fastaStruct(j).Header(1:index(1)-1);
                    [~, index]=ismember(abbrev,phylDistStruct.ids);
                    if any(index)
                        phylDist(j)=phylDistStruct.distMat(index(1),phylDistId);
                    end
                end
            end
            
            %Inf means that it should not be included
            phylDist(phylDist>maxPhylDist)=[];
            
            %Sort based on phylDist
            [~, order]=sort(phylDist);
            
            %Save the first nSequences hits to a temporary FASTA file
            if nSequences<=numel(fastaStruct)
                fastaStruct=fastaStruct(order(1:nSequences));
            else
                fastaStruct=fastaStruct(order);
            end
            
            %Do the alignment if there are more than one sequences,
            %otherwise just save the sequence (or an empty file)
            if numel(fastaStruct)>1
                if (seqIdentity==1 || seqIdentity==0.9 || seqIdentity==0.5)
                    if ~ispc
                        if seqIdentity==1
                            cdhitInp100=tempname;
                            fastawrite(cdhitInp100,fastaStruct);
                            tmpFile=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp100 '" -o "' tmpFile '" -c 1.0 -s 0.8 -n 5 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp100, 'file')
                                delete([cdhitInp100 '*']);
                            end
                        elseif seqIdentity==0.9
                            cdhitInp100=tempname;
                            fastawrite(cdhitInp100,fastaStruct);
                            cdhitInp90=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp100 '" -o "' cdhitInp90 '" -c 1.0 -s 0.8 -n 5 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp100, 'file')
                                delete([cdhitInp100 '*']);
                            end
                            tmpFile=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp90 '" -o "' tmpFile '" -c 0.9 -s 0.8 -n 5 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp90, 'file')
                                delete([cdhitInp90 '*']);
                            end
                        elseif seqIdentity==0.5
                            cdhitInp100=tempname;
                            fastawrite(cdhitInp100,fastaStruct);
                            cdhitInp90=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp100 '" -o "' cdhitInp90 '" -c 1.0 -s 0.8 -n 5 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp100, 'file')
                                delete([cdhitInp100 '*']);
                            end
                            cdhitInp50=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp90 '" -o "' cdhitInp50 '" -c 0.9 -s 0.8 -n 5 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp90, 'file')
                                delete([cdhitInp90 '*']);
                            end
                            tmpFile=tempname;
                            [status, output]=unix(['"' fullfile(ravenPath,'software','cd-hit-v4.6.6',['cd-hit' binEnd]) '" -T "' cores '" -i "' cdhitInp50 '" -o "' tmpFile '" -c 0.5 -s 0.8 -n 3 -M 2000']);
                            if status~=0
                                EM=['Error when performing clustering of ' missingAligned{i} ':\n' output];
                                dispEM(EM);
                            end
                            %Remove the old tempfile
                            if exist(cdhitInp50, 'file')
                                delete([cdhitInp50 '*']);
                            end
                        end
                    else
                        %CD-HIT is not available in Windows so we just
                        %export multifasta file into temporary file;
                        tmpFile=tempname;
                        fastawrite(tmpFile,fastaStruct);
                    end
                else
                    tmpFile=tempname;
                    fastawrite(tmpFile,fastaStruct);
                end
                %Do the alignment for this file
                if ~ispc
                    [status, output]=system(['"' fullfile(ravenPath,'software','mafft-7.305',['mafft' binEnd]) '" --auto --thread "' cores '" "' tmpFile '" > "' fullfile(dataDir,'aligned',[missingAligned{i} '.faw']) '"']);
                else
                    [status, output]=system(['"' fullfile(ravenPath,'software','mafft-7.305','mafft.bat') '" --auto --thread "' cores '" "' tmpFile '" > "' fullfile(dataDir,'aligned',[missingAligned{i} '.faw']) '"']);
                end
                if status~=0
                    EM=['Error when performing alignment of ' missingAligned{i} ':\n' output];
                    dispEM(EM);
                end
                %Remove the old tempfile
                if exist(tmpFile, 'file')
                    delete([tmpFile '*']);
                end
            else
                %If there is only one sequence then it's not possible to do
                %a multiple alignment. Just print the sequence instead. An
                %empty file was written previously so that doesn't have to
                %be dealt with
                if numel(fastaStruct)==1
                    fastawrite(fullfile(dataDir,'aligned',[missingAligned{i} '.faw']),fastaStruct);
                end
            end
            %Move the temporary file to the real one
            movefile(fullfile(dataDir,'aligned',[missingAligned{i} '.faw']),fullfile(dataDir,'aligned',[missingAligned{i} '.fa']),'f');
        end
    end
end

if ~ispc
    if (seqIdentity==1 || seqIdentity==0.9 || seqIdentity==0.5)
        fprintf('Completed clustering and multiple alignment of sequences\n');
    else
        fprintf('Completed multiple alignment of sequences. Protein clustering was not requested or incorrect seqIdentity value was used\n');
    end
else
    if (seqIdentity==1 || seqIdentity==0.9 || seqIdentity==0.5)
        fprintf('Protein clustering was skipped, since CD-HIT is not compatible with Windows');
    end
    fprintf('Completed multiple alignment of sequences\n');
end

%Check if training of Hidden Markov models should be performed
missingHMMs=setdiff(KOModel.rxns,[hmmFiles;outFiles]);
if ~isempty(missingHMMs)
    missingHMMs=missingHMMs(randperm(RandStream.create('mrg32k3a','Seed',cputime()),numel(missingHMMs)));
    
    %Train models for all missing KOs
    for i=1:numel(missingHMMs)
        %This is checked here because it could be that it is created by a
        %parallel process
        if ~exist(fullfile(dataDir,'hmms',[missingHMMs{i} '.hmm']),'file') && ~exist(fullfile(dataDir,'hmms',[missingHMMs{i} '.hmw']),'file')
            %Check that the aligned FASTA file exists. It could be that it
            %is still being worked on by some other instance of the program
            %(the .faw file should then exist). This should not happen on a
            %single computer. It doesn't throw an error, because it should
            %finalize the ones it can
            if ~exist(fullfile(dataDir,'aligned',[missingHMMs{i} '.fa']),'file')
                EM=['The aligned FASTA file for ' missingHMMs{i} ' does not exist'];
                dispEM(EM,false);
                continue;
            end
            
            %If the multi-FASTA file is empty then save an empty aligned
            %file and continue
            s=dir(fullfile(dataDir,'aligned',[missingHMMs{i} '.fa']));
            if s.bytes<=0
                fid=fopen(fullfile(dataDir,'hmms',[missingHMMs{i} '.hmm']),'w');
                fclose(fid);
                continue;
            end
            %Create a temporary file to indicate that it is working on the
            %KO. This is because hmmbuild cannot overwrite existing files
            fid=fopen(fullfile(dataDir,'hmms',[missingHMMs{i} '.hmw']),'w');
            fclose(fid);
            
            %Create HMM
            [status, output]=system(['"' fullfile(ravenPath,'software','hmmer-3.1b2',['hmmbuild' binEnd]) '" --cpu "' cores '" "' fullfile(dataDir,'hmms',[missingHMMs{i} '.hmm']) '" "' fullfile(dataDir,'aligned',[missingHMMs{i} '.fa']) '"']);
            if status~=0
                EM=['Error when training HMM for ' missingHMMs{i} ':\n' output];
                dispEM(EM);
            end
            
            %Delete the temporary file
            delete(fullfile(dataDir,'hmms',[missingHMMs{i} '.hmw']));
        end
    end
end
fprintf('Completed generation of HMMs\n');

%Check which new .out files that should be generated
%Check if training of Hidden Markov models should be performed
missingOUT=setdiff(KOModel.rxns,outFiles);
if ~isempty(missingOUT)
    missingOUT=missingOUT(randperm(RandStream.create('mrg32k3a','Seed',cputime()),numel(missingOUT)));
    for i=1:numel(missingOUT)
        %This is checked here because it could be that it is created by a
        %parallel process
        if ~exist(fullfile(outDir,[missingOUT{i} '.out']),'file')
            %Check that the HMM file exists. It should do so since %we are
            %saving empty files as well. Print a warning and continue if
            %not
            if ~exist(fullfile(dataDir,'hmms',[missingOUT{i} '.hmm']),'file')
                EM=['The HMM file for ' missingOUT{i} ' does not exist'];
                dispEM(EM,false);
                continue;
            end
            
            %Save an empty file to prevent several threads working on the
            %same file
            fid=fopen(fullfile(outDir,[missingOUT{i} '.out']),'w');
            fclose(fid);
            
            %If the HMM file is empty then save an out file and continue
            s=dir(fullfile(dataDir,'hmms',[missingOUT{i} '.hmm']));
            if s.bytes<=0
                continue;
            end
            
            %Check each gene in the input file against this model
            [status, output]=system(['"' fullfile(ravenPath,'software','hmmer-3.1b2',['hmmsearch' binEnd]) '" --cpu "' cores '" "' fullfile(dataDir,'hmms',[missingOUT{i} '.hmm']) '" "' fastaFile '"']);
            if status~=0
                EM=['Error when querying HMM for ' missingOUT{i} ':\n' output];
                dispEM(EM);
            end
            
            %Save the output to a file
            fid=fopen(fullfile(outDir,[missingOUT{i} '.out']),'w');
            fwrite(fid,output);
            fclose(fid);
        end
    end
end
fprintf('Completed matching to HMMs\n');

%***Begin retrieving the output and putting together the resulting model

%Retrieve matched genes from the HMMs
koGeneMat=zeros(numel(KOModel.rxns),3000); %Make room for 3000 genes
genes=cell(3000,1);
%Store the best score for a gene in a hash list (since it will be searching
%many times)
hTable = java.util.Hashtable;

geneCounter=0;
for i=1:numel(KOModel.rxns)
    if exist(fullfile(outDir,[KOModel.rxns{i} '.out']), 'file')
        fid=fopen(fullfile(outDir,[KOModel.rxns{i} '.out']),'r');
        beginMatches=false;
        while 1
            %Get the next line
            tline = fgetl(fid);
            
            %Abort at end of file
            if ~ischar(tline)
                break;
            end
            
            if and(beginMatches,strcmp(tline,'  ------ inclusion threshold ------'))
                break;
            end
            
            if beginMatches==false
                %This is how the listing of matches begins
                if any(strfind(tline,'E-value '))
                    %Read one more line that is only padding
                    tline = fgetl(fid);
                    beginMatches=true;
                end
            else
                %If matches should be read
                if ~strcmp(tline,'   [No hits detected that satisfy reporting thresholds]') && ~isempty(tline)
                    elements=regexp(tline,' ','split');
                    elements=elements(cellfun(@any,elements));
                    
                    %Check if the match is below the treshhold
                    score=str2num(elements{1});
                    gene=elements{9};
                    if score<=cutOff
                        %If the score is exactly 0, change it to a very
                        %small value to avoid NaN
                        if score==0
                            score=10^-250;
                        end
                        %Check if the gene is added already and, is so, get
                        %the best score for it
                        I=hTable.get(gene);
                        if any(I)
                            koGeneMat(i,I)=score;
                        else
                            geneCounter=geneCounter+1;
                            %The gene was not present yet so add it
                            hTable.put(gene,geneCounter);
                            genes{geneCounter}=gene;
                            koGeneMat(i,geneCounter)=score;
                        end
                    end
                else
                    break;
                end
            end
        end
        fclose(fid);
    end
end
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

%Find and delete all reactions without any genes. This also removes genes
%that are not used (which could happen because minScoreRatioG and
%minScoreRatioKO)
I=sum(model.rxnGeneMat,2)==0;
model=removeReactions(model,I,true,true);

%Add the gene associations as 'or'
for i=1:numel(model.rxns)
    %Find the involved genes
    I=find(model.rxnGeneMat(i,:));
    model.grRules{i}=['(' model.genes{I(1)}];
    for j=2:numel(I)
        model.grRules{i}=[model.grRules{i} ' or ' model.genes{I(j)}];
    end
    model.grRules{i}=[model.grRules{i} ')'];
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
end

%Supporter function to list the files in a directory and return them as a
%cell array
function files=listFiles(directory)
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
