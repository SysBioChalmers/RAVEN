function [draftModel, hitGenes]=getModelFromHomology(models,blastStructure,...
    getModelFor,preferredOrder,strictness,onlyGenesInModels,maxE,...
    minLen,minIde,mapNewGenesToOld)
% getModelFromHomology
%   Constructs a new model from a set of existing models and gene homology
%   information.
%
%   models            a cell array of model structures to build the model
%                     from. These models must be sorted by importance in
%                     decreasing order
%   blastStructure    a blastStructure as produced by getBlast or
%                     getBlastFromExcel
%   getModelFor       a three-four letter abbreviation of the organism to
%                     build a model for. Must have BLASTP hits in both
%                     directions to the organisms in 'models'
%   preferredOrder    the order in which reactions should be added from the
%                     models. If not supplied, reactions will be included
%                     from all models, otherwise one gene will only result
%                     in reactions from one model (opt, default {})
%   strictness        integer that specifies which reactions should be
%                     included:
%                     1: Map new genes to old for all pairs, which have
%                     acceptable BLASTP results in both directions
%                     2: Map new genes to old for all pairs, which have
%                     acceptable BLASTP results in correspondent direction
%                     (mapping can be done in the opposite direction, see
%                     mapNewGenesToOld below)
%                     3: Check all BLASTP results and retain only the best
%                     results by E-value for all gene pairs in each
%                     direction separately. Then map new genes to old for
%                     all pairs, which have acceptable BLASTP results in
%                     both directions (opt, default 1).
%   onlyGenesInModels consider BLASTP results only for genes that exist in
%                     the models. This tends to import a larger fraction
%                     from the existing models but may give less reliable
%                     results. Has effect only if strictness=3 (opt,
%                     default false)
%   maxE              only look at genes with E-values <= this value (opt,
%                     default 10^-30)
%   minLen            only look at genes with alignment length >= this
%                     value (opt, default 200)
%   minIde            only look at genes with identity >= this value
%                     (opt, default 40 (%))
%   mapNewGenesToOld  determines how to match genes if not looking at only
%                     1-1 orthologs. Either map the new genes to the old or
%                     old genes to new. The default is to map the new genes
%                     (opt, default true)
%
%   draftModel        a model structure for the new organism
%   hitGenes          collect the old and new genes
%
%   The models in the 'models' structure should have named the metabolites
%   in the same manner, have their reversible reactions in the same
%   direction (run sortModel), and use the same compartment names. To avoid
%   keeping unneccesary old genes, the models should not have
%   'or'-relations in their grRules (use expandModel).
%
%   The resulting draft model contains only reactions associated with
%   orthologous genes. The old (original) genes involved in 'and'
%   relations in grRules without any orthologs are still included in
%   the draft model as OLD_MODELID_geneName.
%
%   NOTE: "to" and "from" means relative to the new organism
%
%   Usage: draftModel=getModelFromHomology(models,blastStructure,...
%    getModelFor,preferredOrder,strictness,onlyGenesInModels,maxE,...
%    minLen,minIde,mapNewGenesToOld)

hitGenes.oldGenes = [];  % collect the old genes from the template model (organism)
hitGenes.newGenes = [];  % collect the new genes of the draft model (target organism)

getModelFor=char(getModelFor);

if nargin<4
    preferredOrder=[];
else
    preferredOrder=convertCharArray(preferredOrder);
    [row,col]=size(preferredOrder);
    if col>row
        preferredOrder=transpose(preferredOrder);
    end
end
if nargin<5
    strictness=1;
end
if nargin<6
    onlyGenesInModels=false;
end
if nargin<7
    maxE=10^-30;
end
if nargin<8
    minLen=200;
end
if nargin<9
    minIde=40;
end
if nargin<10
    mapNewGenesToOld=true;
end

if isfield(models,'S')
    models={models};
end
    
%Check that all the information is in the blast structure
modelNames=cell(numel(models),1);
for i=1:numel(models)
    modelNames{i}=models{i}.id;
    %Gene short names and geneMiriams are often different between species,
    %safer not to include them
    if isfield(models{i},'geneShortNames')
        models{i}=rmfield(models{i},'geneShortNames');
    end
    if isfield(models{i},'geneMiriams')
        models{i}=rmfield(models{i},'geneMiriams');
    end
end

%Assume for now that all information is there and that it's correct. This
%is important to fix since no further checks are being made!

%Check whether provided fasta files use the same gene identifiers as
%provided template models
for i=1:numel(blastStructure)
    if ~strcmp(blastStructure(i).fromId,getModelFor)
        j=strcmpi(blastStructure(i).fromId,modelNames);
        if j==0
            error(['While the blastStructure contains sequences from '...
                'organismID "%s" (as\nprovided in getBlast), none of '...
                'template models have this id (as model.id)'],...
                string(blastStructure(i).fromId));
        end
        k=sum(ismember(blastStructure(i).fromGenes,models{j}.genes));
        if k<(numel(models{j}.genes)*0.05)
            error(['Less than 5%% of the genes in the template model '...
                'with model.id "%s"\ncan be found in the blastStructure. '...
                'Ensure that the protein FASTA\nused in getBlast and '...
                'the template model used in getModelFromHomology\nuse '...
                'the same style of gene identifiers'],models{j}.id)
        end
    end
end

%Standardize grRules of template models
for i=1:length(models)
    fprintf('\nStandardizing grRules of template model with ID "%s" ...',models{i}.id);
    [models{i}.grRules,models{i}.rxnGeneMat]=standardizeGrRules(models{i},false);
end
fprintf(' done\n');

%Remove all gene matches that are below the cutoffs
for i=1:numel(blastStructure)
    indexes=blastStructure(i).evalue<maxE & blastStructure(i).aligLen>=minLen & blastStructure(i).identity>=minIde; %Do it in this direction to lose NaNs
    blastStructure(i).fromGenes(~indexes)=[];
    blastStructure(i).toGenes(~indexes)=[];
    blastStructure(i).evalue(~indexes)=[];
    blastStructure(i).identity(~indexes)=[];
    blastStructure(i).aligLen(~indexes)=[];
    blastStructure(i).bitscore(~indexes)=[];
    blastStructure(i).ppos(~indexes)=[];
end

%Remove all reactions from the models that have no genes encoding for them.
%Also remove all genes that encode for no reactions. There shouldn't be any
%but there might be mistakes
for i=1:numel(models)
    [hasGenes, ~]=find(models{i}.rxnGeneMat);
    hasNoGenes=1:numel(models{i}.rxns);
    hasNoGenes(hasGenes)=[];
    models{i}=removeReactions(models{i},hasNoGenes,true,true);
end

%Create a structure that contains all genes used in the blasts in any
%direction for each of the models in 'models' and for the new organism. The
%first cell is for the new organism and then according to the preferred
%order. If no such order is supplied, then according to the order in
%'models'
allGenes=cell(numel(models)+1,1);
if isempty(preferredOrder)
    useOrder=modelNames;
else
    useOrder=preferredOrder;
end

%Get the corresponding indexes for those models in the 'models' structure
useOrderIndexes=zeros(numel(models),1);
for i=1:numel(models)
    [~, index]=ismember(models{i}.id,useOrder);
    useOrderIndexes(index)=i;
end

%Remove all genes from the blast structure that have no genes in the models
if onlyGenesInModels==true
    modelGenes={};
    for i=1:numel(models)
        modelGenes=[modelGenes;models{i}.genes];
    end
    for i=1:numel(blastStructure)
        %Check to see if it should match the toId or fromId
        if strcmpi(blastStructure(i).fromId,getModelFor)
            I=ismember(blastStructure(i).toGenes,modelGenes);
        else
            I=ismember(blastStructure(i).fromGenes,modelGenes);
        end
        blastStructure(i).fromGenes(~I)=[];
        blastStructure(i).toGenes(~I)=[];
        blastStructure(i).evalue(~I)=[];
        blastStructure(i).identity(~I)=[];
        blastStructure(i).aligLen(~I)=[];
        blastStructure(i).bitscore(~I)=[];
        blastStructure(i).ppos(~I)=[];
        
        %Check that no matching in blastStructure is empty. This happens if
        %no genes in the models are present in the corresponding sheet
        if isempty(blastStructure(i).fromGenes)
            EM=['No genes in matching from ' blastStructure(i).fromId ' to ' blastStructure(i).toId ' are present in the corresponding model'];
            dispEM(EM);
        end
    end
end

%If only best orthologs are to be used then all other measurements are
%deleted from the blastStructure. All code after this stays the same. This
%means that preferred order can still matter. The best ortholog scoring is
%based only on the E-value
if strictness==3
    for i=1:numel(blastStructure)
        keep=false(numel(blastStructure(i).toGenes),1);
        [allFromGenes, ~, I]=unique(blastStructure(i).fromGenes);
        
        %It would be nice to get rid of this loop
        for j=1:numel(allFromGenes)
            allMatches=find(I==j);
            bestMatches=allMatches(blastStructure(i).evalue(allMatches)==min(blastStructure(i).evalue(allMatches)));
            
            %Keep the best matches
            keep(bestMatches)=true;
        end
        
        %Delete all matches that were not best matches
        blastStructure(i).fromGenes(~keep)=[];
        blastStructure(i).toGenes(~keep)=[];
        blastStructure(i).evalue(~keep)=[];
        blastStructure(i).identity(~keep)=[];
        blastStructure(i).aligLen(~keep)=[];
        blastStructure(i).bitscore(~keep)=[];
        blastStructure(i).ppos(~keep)=[];
    end
end

useOrder=[{getModelFor};useOrder];

for i=1:numel(blastStructure)
    [~, toIndex]=ismember(blastStructure(i).toId,useOrder);
    [~, fromIndex]=ismember(blastStructure(i).fromId,useOrder);
    
    %Add all genes to the corresponding list in allGenes
    allGenes{toIndex}=[allGenes{toIndex};blastStructure(i).toGenes];
    allGenes{fromIndex}=[allGenes{fromIndex};blastStructure(i).fromGenes];
end

%Keep only the unique gene names
maxOtherGeneNr=0; %Determines the dimension of the connectivity matrixes
for i=1:numel(allGenes)
    allGenes{i}=unique(allGenes{i});
    if i>1
        if numel(allGenes{i})>maxOtherGeneNr
            maxOtherGeneNr=numel(allGenes{i});
        end
    end
end

%Generate a cell array of matrixes that describes how the genes in the new
%organism map to the models. Each cell matches to the corresponding model
%in useOrder (starting at 2 of course). First dimension is gene in new
%organism, second which gene it is in the other organism. The second matrix
%describes how they map back.

%As it is now, a significant match is indicated by a 1. This could be
%expanded to contain some kind of significance level. The first dimension
%is still the genes in the new model.

allTo=cell(numel(useOrder)-1,1);
allFrom=cell(numel(useOrder)-1,1);

for i=1:numel(useOrder)-1
    allTo{i}=sparse(numel(allGenes{1}),numel(allGenes{i+1}));
    allFrom{i}=sparse(numel(allGenes{1}),numel(allGenes{i+1}));
end

%Fill the matches to other species
for i=1:numel(blastStructure)
    if strcmp(blastStructure(i).toId,getModelFor)
        %This was 'to' the new organism. They should all match so no checks
        %are being made
        [~, a]=ismember(blastStructure(i).toGenes,allGenes{1});
        [~, fromModel]=ismember(blastStructure(i).fromId,useOrder);
        [~, b]=ismember(blastStructure(i).fromGenes,allGenes{fromModel});
        idx = sub2ind(size(allTo{fromModel-1}), a, b);
        allTo{fromModel-1}(idx)=1;
    else
        %This was 'from' the new organism
        [~, a]=ismember(blastStructure(i).fromGenes,allGenes{1});
        [~, toModel]=ismember(blastStructure(i).toId,useOrder);
        [~, b]=ismember(blastStructure(i).toGenes,allGenes{toModel});
        idx = sub2ind(size(allFrom{toModel-1}), a, b);
        allFrom{toModel-1}(idx)=1;
    end
end

%Now we have all the gene matches in a convenient way. For all the genes in
%the new organism get the genes that should be included from other
%organisms. If all genes should be included this simply means keep the
%allFrom matrix as it is. If only orthologs which could be mapped in both
%BLAST directions are to be included then only those elements are kept.

finalMappings=cell(numel(useOrder)-1,1);
if strictness==1 || strictness==3
    for j=1:numel(allFrom)
        finalMappings{j}=allTo{j}~=0 & allFrom{j}~=0;
    end
else
    if mapNewGenesToOld==true
        finalMappings=allFrom;
    else
        finalMappings=allTo;
    end
end

%Remove all genes from the mapping that are not in the models. This doesn't
%do much if only genes in the models were used for the original mapping.
%Also simplify the finalMapping and allGenes structures so that they only
%contain mappings that exist
usedNewGenes=false(numel(allGenes{1}),1);

for i=1:numel(allGenes)-1
    %First remove mappings for those genes that are not in the model
    if onlyGenesInModels==false
        a=ismember(allGenes{i+1},models{useOrderIndexes(i)}.genes);
        finalMappings{i}(:,~a)=false;
    end
    
    %Then remove unused ones and simplify
    [a, b]=find(finalMappings{i});
    usedGenes=false(numel(allGenes{i+1}),1);
    usedNewGenes(a)=true;
    usedGenes(b)=true;
    finalMappings{i}=finalMappings{i}(:,usedGenes);
    allGenes{i+1}=allGenes{i+1}(usedGenes);
end

%Remove all new genes that have not been mapped to anything
allGenes{1}=allGenes{1}(usedNewGenes);
for i=1:numel(finalMappings)
    finalMappings{i}=finalMappings{i}(usedNewGenes,:);
end

%Now is it time to choose which reactions should be included from which
%models. If there is a preferred order specified then each gene can only
%result in reactions from one model, otherwise they should all be included

%Start by simplifying the models by removing genes/reactions that are not
%used. This is where it gets weird with complexes, especially "or"
%complexes. In this step only reactions which are encoded by one single
%gene, or where all genes should be deleted, are deleted. The info on the
%full complex is still present in the grRules

for i=1:numel(models)
    a=ismember(models{useOrderIndexes(i)}.genes,allGenes{i+1});
    
    %Remove reactions that are not associated to any of the genes in
    %allGenes, thereby also keeping complexes where only for one of the
    %genes was matched
    [rxnsToKeep,~] = find(models{useOrderIndexes(i)}.rxnGeneMat(:,a));
    rxnsToRemove = repmat(1,numel(models{useOrderIndexes(i)}.rxns),1);
    rxnsToRemove(rxnsToKeep) = 0;
    rxnsToRemove = find(rxnsToRemove);
    models{useOrderIndexes(i)}=removeReactions(models{useOrderIndexes(i)},rxnsToRemove,true,true,true);
end

%Since mergeModels function will be used in the end, the models are
%simplified further by deleting genes/reactions in the order specified by
%preferredOrder. This means that the last model will only contain reactions
%for genes that mapped only to that model

allUsedGenes=false(numel(allGenes{1}),1);

if ~isempty(preferredOrder) && numel(models)>1
    [usedGenes, ~]=find(finalMappings{1}); %All that are used in the first model in preferredOrder
    allUsedGenes(usedGenes)=true;
    for i=2:numel(finalMappings)
        [usedGenes, ~]=find(finalMappings{i});
        usedGenes=unique(usedGenes);
        a=ismember(usedGenes,find(allUsedGenes));
        
        [~, genesToDelete]=find(finalMappings{i}(usedGenes(a),:)); %IMPORTANT! IS it really correct to remove all genes that map?
        genesToDelete=unique(genesToDelete); %Maybe not needed, but for clarity if nothing else
        
        %Remove all the genes that were already found and add the other
        %ones to allUsedGenes
        models{useOrderIndexes(i)}=removeGenes(models{useOrderIndexes(i)},allGenes{i+1}(genesToDelete),true,true,false);
        allUsedGenes(usedGenes)=true;
        
        %Remove the deleted genes from finalMappings and allGenes.
        finalMappings{i}(:,genesToDelete)=[];
        allGenes{i+1}(genesToDelete)=[];
    end
end

%Now loop through the models and update the gene associations. Genes not
%belonging to the new organism will be renamed as 'OLD_MODELID_gene'
for i=1:numel(models)
    %Find all the new genes that should be used for this model
    [newGenes, oldGenes]=find(finalMappings{i});
    
    %Create a new gene list with the genes from the new organism and those
    %genes that could not be removed
    replaceableGenes=allGenes{i+1}(unique(oldGenes));
    nonReplaceableGenes=setdiff(models{useOrderIndexes(i)}.genes,replaceableGenes);
    fullGeneList=[allGenes{1}(unique(newGenes));nonReplaceableGenes];
    
    %Just to save some indexing later. This is the LAST index of
    %replaceable ones
    nonRepStartIndex=numel(unique(newGenes));
    
    %Construct a new rxnGeneMat
    newRxnGeneMat=sparse(numel(models{useOrderIndexes(i)}.rxns),numel(fullGeneList));
    
    %Now update the rxnGeneMat. This is a little tricky and could
    %probably be done in a more efficient way, but I just loop through the
    %reactions and add them one by one
    for j=1:numel(models{useOrderIndexes(i)}.rxns)
        %Get the old genes encoding for this reaction
        [~, oldGeneIds]=find(models{useOrderIndexes(i)}.rxnGeneMat(j,:));
        
        %Update the matrix for each gene. This includes replacing one gene
        %with several new ones if there were several matches
        for k=1:numel(oldGeneIds)
            %Match the gene to one in the gene list. This is done as a text
            %match. Could probably be done better, but I'm a little lost in
            %the indexing
            
            geneName=models{useOrderIndexes(i)}.genes(oldGeneIds(k));
            
            %First search in the mappable genes
            mapIndex=find(ismember(allGenes{i+1},geneName));
            
            if ~isempty(mapIndex)
                % add the old genes
                hitGenes.oldGenes = [hitGenes.oldGenes, {geneName}];
                
                %Get the new genes for that gene
                a=find(finalMappings{i}(:,mapIndex));
                
                %Find the positions of these genes in the final gene list
                [~, b]=ismember(allGenes{1}(a),fullGeneList);
                
                %Update the matrix
                newRxnGeneMat(j,b)=1;
                
                %Update the grRules string. This is tricky, but I hope that
                %it's ok to replace the old gene name with the new one and
                %add ') or (' if there were several matches. Be sure of
                %this!
                repString=fullGeneList{b(1)};
                if numel(b)>1
                    for l=2:numel(b)
                        repString=[repString ' or ' fullGeneList{b(l)}];
                    end
                    repString=['(' repString ')'];
                end
                
                % add the new matched genes
                hitGenes.newGenes = [hitGenes.newGenes, {repString}];
                
                %Use regexprep instead of strrep to prevent partial matches
                models{useOrderIndexes(i)}.grRules{j}=regexprep(models{useOrderIndexes(i)}.grRules{j},['(^|\s|\()' geneName{1} '($|\s|\))'],['$1' repString '$2']);
            else
                %Then search in the non-replaceable genes. There could only
                %be one match here
                index=find(ismember(nonReplaceableGenes,geneName));
                
                %Update the matrix
                newRxnGeneMat(j,nonRepStartIndex+index)=1;
                
                models{useOrderIndexes(i)}.grRules{j}=strrep(models{useOrderIndexes(i)}.grRules{j},geneName{1},strcat('OLD_',models{useOrderIndexes(i)}.id,'_',geneName{1}));
            end
        end
    end
    
    %Add the new list of genes
    models{useOrderIndexes(i)}.rxnGeneMat=newRxnGeneMat;
    if ~isempty(nonReplaceableGenes)
        models{useOrderIndexes(i)}.genes=[allGenes{1}(unique(newGenes));strcat('OLD_',models{useOrderIndexes(i)}.id,'_',nonReplaceableGenes)];
    else
        models{useOrderIndexes(i)}.genes=allGenes{1}(unique(newGenes));
    end
    if isfield(models{useOrderIndexes(i)},'geneComps')
        geneComps=models{useOrderIndexes(i)}.geneComps(1);
        models{useOrderIndexes(i)}.geneComps=zeros(numel(models{useOrderIndexes(i)}.genes),1);
        %Assume that all genes are in the same compartment, and this
        %compartment is specified for the first gene
        models{useOrderIndexes(i)}.geneComps(:)=geneComps;
    end
end

%Now merge the models. All information should be correct except for 'or'
%complexes
draftModel=mergeModels(models,'metNames');

% %Remove unnecessary OLD_ genes, that were added with OR relationships
regexStr=['OLD_(', strjoin(modelNames(:),'|'),')_\S+'];
draftModel.grRules=regexprep(draftModel.grRules,[' or ' regexStr],'');
draftModel.grRules=regexprep(draftModel.grRules,[regexStr ' or '],'');

%Change name of the resulting model
draftModel.id=getModelFor;
name='Generated by getModelFromHomology using ';
for i=1:numel(models)
    if i<numel(models)
        name=[name models{i}.id ', '];
    else
        name=[name models{i}.id];
    end
end
draftModel.name=name;
draftModel.rxnNotes=cell(length(draftModel.rxns),1);
draftModel.rxnNotes(:)={'Included by getModelFromHomology'};
draftModel.rxnConfidenceScores=NaN(length(draftModel.rxns),1);
draftModel.rxnConfidenceScores(:)=2;
draftModel=deleteUnusedGenes(draftModel,0);
%Standardize grRules and notify if problematic grRules are found
[draftModel.grRules,draftModel.rxnGeneMat]=standardizeGrRules(draftModel,false);
draftModel=deleteUnusedGenes(draftModel,false);
end
