function [model,KOModel]=getModelFromKEGG(keggPath,keepSpontaneous,...
    keepUndefinedStoich,keepIncomplete,keepGeneral)
% getModelFromKEGG
%   Retrieves information stored in KEGG database and generates a model
%
%   Input:
%   keggPath            if keggGenes.mat, keggMets.mat, keggPhylDist.mat or
%                       keggRxns.mat are not in the RAVEN/external/kegg
%                       directory, this function will attempt to read data
%                       from a local FTP dump of the KEGG database.
%                       keggPath is the path to the root of this database
%                       (opt, default 'RAVEN/external/kegg'). If
%                       keggModel.mat is present in the same directory, the
%                       function reads the data from this file and ignores
%                       keggGenes.mat, keggMets.mat and keggRxns.mat
%   keepSpontaneous     include reactions labeled as "spontaneous" (opt,
%                       default true)
%   keepUndefinedStoich include reactions in the form n A <=> n+1 A. These
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
%
%   Output:
%   model               a model structure generated from the database. All
%                       reactions and the metabolites used in them will be
%                       added
%   KOModel             a model structure representing the KEGG Orthology
%                       ids and their associated genes. The KO ids are
%                       saved as reactions
%
%   NOTE: The model output from getModelFromKEGG can be used as template
%   for fillGaps. In that case, ensure that the genes and rxnGeneMat fields
%   are removed before parsing: model=rmfield(model,'genes'), etc.
%
%   Usage: [model,KOModel]=getModelFromKEGG(keggPath,keepSpontaneous,...
%    keepUndefinedStoich,keepIncomplete,keepGeneral)

if nargin<1
    keggPath='RAVEN/external/kegg';
else
    keggPath=char(keggPath);
end
if nargin<2
    keepSpontaneous=true;
end
if nargin<3
    keepUndefinedStoich=true;
end
if nargin<4
    keepIncomplete=true;
end
if nargin<5
    keepGeneral=false;
end

ravenPath=findRAVENroot();
modelFile=fullfile(ravenPath,'external','kegg','keggModel.mat');
if exist(modelFile, 'file') && isNewestFile(ravenPath)
    fprintf(['Importing the global KEGG model from ' strrep(modelFile,'\','/') '... ']);
    load(modelFile);
    fprintf('COMPLETE\n');
else
    fprintf(['NOTE: The file ' strrep(modelFile,'\','/') ' does not exist or is out-of-date and therefore will be (re)generated\n']);
    %First get all reactions
    [model,isSpontaneous,isUndefinedStoich,isIncomplete,isGeneral]=getRxnsFromKEGG(keggPath);
    
    %Get the KO ids that are associated with any of the reactions. They
    %will be used later on to create a rxn-gene matrix
    KOs=cell(numel(model.rxns)*2,1);
    %Make room for two KO ids per reaction
    
    addToIndex=1;
    %Add all KO ids that are used
    for i=1:numel(model.rxns)
        if isstruct(model.rxnMiriams{i})
            for j=1:numel(model.rxnMiriams{i}.name)
                if strcmpi('kegg.orthology',model.rxnMiriams{i}.name{j})
                    %Add the KO id
                    KOs(addToIndex)=model.rxnMiriams{i}.value(j);
                    addToIndex=addToIndex+1;
                end
            end
        end
    end
    
    KOs=KOs(1:addToIndex-1);
    KOs=unique(KOs);
    
    %Get all genes from any organism in KEGG that is associated with any of
    %the KOs
    KOModel=getGenesFromKEGG(keggPath,KOs);
    
    fprintf('Pruning the global KEGG model from the partially annotated, lumped KEGG Orthology entries... ')
    model.genes=KOModel.genes;
    
    %It can be that there are KOs from the reactions that have no database
    %entry. These are (as far as I have seen) lumped versions of other KOs
    %and should be removed
    KOsToRemove=setdiff(KOs, KOModel.rxns);
    
    %Loop through all reactions and delete the KOs that were not found
    for i=1:numel(model.rxns)
        if isstruct(model.rxnMiriams{i})
            for j=1:numel(model.rxnMiriams{i}.name)
                toDel=[];
                if strcmp(model.rxnMiriams{i}.name{j},'kegg.orthology')
                    if ismember(model.rxnMiriams{i}.value{j},KOsToRemove)
                        toDel=[toDel;j];
                    end
                end
            end
            %Delete the KOs
            if any(toDel)
                %If all posts are deleted, then delete the whole structure
                if numel(toDel)==j
                    model.rxnMiriams{i}=[];
                else
                    model.rxnMiriams{i}.name(toDel)=[];
                    model.rxnMiriams{i}.value(toDel)=[];
                end
            end
        end
    end
    fprintf('COMPLETE\n');
    
    fprintf('Constructing the rxnGeneMat for the global KEGG model...   0%% complete');
    %Create the rxnGeneMat for the reactions. This is simply done by
    %merging the gene associations for all the involved KOs
    r=zeros(10000000,1);
    %Store the positions since it's slow to write to a sparse array in a
    %loop
    c=zeros(10000000,1);
    counter=1;
    for i=1:numel(model.rxns)
        if isstruct(model.rxnMiriams{i})
            I=strncmp('kegg.orthology',model.rxnMiriams{i}.name,18);
            if any(I)
                [J, K]=ismember(model.rxnMiriams{i}.value(I),KOModel.rxns);
                %Find all gene indexes that correspond to any of these KOs
                [~, L]=find(KOModel.rxnGeneMat(K(J),:));
                if any(L)
                    %Allocate room for more elements if needed
                    if counter+numel(L)-1>=numel(r)
                        r=[r;zeros(numel(r),1)];
                        c=[c;zeros(numel(c),1)];
                    end
                    r(counter:counter+numel(L)-1)=ones(numel(L),1)*i;
                    c(counter:counter+numel(L)-1)=L(:);
                    counter=counter+numel(L);
                end
            end
        end
        if rem(i-1,100) == 0
            progress=pad(num2str(floor(i/numel(model.rxns)*100)),3,'left');
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
        end
    end
    
    model.rxnGeneMat=sparse(r(1:counter-1),c(1:counter-1),ones(counter-1,1));
    if size(model.rxnGeneMat,1)~=numel(model.rxns) || size(model.rxnGeneMat,2)~=numel(KOModel.genes)
        model.rxnGeneMat(numel(model.rxns),numel(KOModel.genes))=0;
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');
    
    %Then get all metabolites
    metModel=getMetsFromKEGG(keggPath);
    
    fprintf('Finalizing the global KEGG model... ');
    %Add information about all metabolites to the model
    [a, b]=ismember(model.mets,metModel.mets);
    a=find(a);
    b=b(a);
    
    if ~isfield(model,'metNames')
        model.metNames=cell(numel(model.mets),1);
        model.metNames(:)={''};
    end
    model.metNames(a)=metModel.metNames(b);
    
    if ~isfield(model,'metFormulas')
        model.metFormulas=cell(numel(model.mets),1);
        model.metFormulas(:)={''};
    end
    model.metFormulas(a)=metModel.metFormulas(b);
    
    if ~isfield(model,'inchis')
        model.inchis=cell(numel(model.mets),1);
        model.inchis(:)={''};
    end
    model.inchis(a)=metModel.inchis(b);
    
    if ~isfield(model,'metMiriams')
        model.metMiriams=cell(numel(model.mets),1);
    end
    model.metMiriams(a)=metModel.metMiriams(b);
    
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
    
    %Put all metabolites in one compartment called 's' (for system). This
    %is done just to be more compatible with the rest of the code
    model.comps={'s'};
    model.compNames={'System'};
    model.metComps=ones(numel(model.mets),1);
    
    %If reactions with undefined stoichiometry are kept, then the
    %corresponding metabolites will have ids such as "(n+1) C00001" and
    %their names will be empty. These ids are not valid SBML identifiers
    %and are therefore replaced with "undefined1, undefined2...". The
    %former ids are stored as the new names
    I=find(cellfun(@any,strfind(model.mets,'n')) | cellfun(@any,strfind(model.mets,'m')));
    model.metNames(I)=model.mets(I);
    repNums=1:numel(I);
    repIDs=strcat('undefined_',cellfun(@num2str,num2cell(repNums(:)),'UniformOutput',false));
    model.mets(I)=repIDs;
    
    %It could also be that the metabolite names are empty for some other
    %reason. In that case, use the ID instead
    I=cellfun(@isempty,model.metNames);
    model.metNames(I)=model.mets(I);

    %Deafult LB and UB
    model.annotation.defaultLB = -1000;
    model.annotation.defaultUB = 1000;
    
    %Save the model structure
    save(modelFile,'model','KOModel','isGeneral','isIncomplete','isUndefinedStoich','isSpontaneous');
    fprintf('COMPLETE\n');
end

%Delete reactions which are labeled as "incomplete", "erroneous",
%"unclear", "general reaction" or having undefined stoichiometry (depending
%on settings)
if keepSpontaneous==false
    model=removeReactions(model,intersect(isSpontaneous,model.rxns),true,true);
end
if keepUndefinedStoich==false
    model=removeReactions(model,intersect(isUndefinedStoich,model.rxns),true,true);
end
if keepIncomplete==false
    model=removeReactions(model,intersect(isIncomplete,model.rxns),true,true);
end
if keepGeneral==false
    model=removeReactions(model,intersect(isGeneral,model.rxns),true,true);
end

end

function output = isNewestFile(ravenPath)
%The ad hoc function, which checks whether keggModel.mat is the more
%recently modified than keggRxns.mat, keggGenes.mat and keggRxns.mat
modelFile=fullfile(ravenPath,'external','kegg','keggModel.mat');
rxnsFile=fullfile(ravenPath,'external','kegg','keggRxns.mat');
genesFile=fullfile(ravenPath,'external','kegg','keggGenes.mat');
metsFile=fullfile(ravenPath,'external','kegg','keggMets.mat');
if (getFileTime(modelFile)>getFileTime(rxnsFile))&&...
        (getFileTime(modelFile)>getFileTime(genesFile))&&...
        (getFileTime(modelFile)>getFileTime(metsFile))
    output=1;
else
    output=0;
end
end

function modTime = getFileTime(fileName)
%Gets a last modification time for a particular file in datenum format that
%the numbers could be easily compared for different files
listing = dir(fileName);
assert(numel(listing) == 1, 'No such file: %s', fileName);
modTime = listing.datenum;
format long;
end
