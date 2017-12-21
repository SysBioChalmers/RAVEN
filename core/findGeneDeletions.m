function [genes, fluxes, originalGenes, details]=findGeneDeletions(model,testType,analysisType,refModel,oeFactor)
% findGeneDeletions
%   Deletes genes, optimizes the model, and keeps track of the resulting
%   fluxes. This is used for identifying gene deletion targets.
%
%   model           a model structure
%   testType        single/double gene deletions/over expressions. Over
%                   expression only available if using MOMA
%                   'sgd'   single gene deletion
%                   'dgd'   double gene deletion
%                   'sgo'   single gene over expression
%                   'dgo'   double gene over expression
%   analysisType    determines whether to use FBA ('fba') or MOMA ('moma')
%                   in the optimization
%   refModel        MOMA works by fitting the flux distributions of two
%                   models to be as similar as possible. The most common
%                   application is where you have a reference model where
%                   some of the fluxes are constrained from experimental
%                   data. This model is required when using MOMA
%   oeFactor        a factor by which the fluxes should be increased if a
%                   gene is overexpressed (opt, default 10)
%
%   genes           a matrix with the genes that were deleted in each
%                   optimization (the gene indexes in originalGenes). Each
%                   row corresponds to a column in fluxes
%   fluxes          a matrix with the resulting fluxes. Double deletions
%                   that result in an unsolvable problem have all zero
%                   flux. Single deletions that result in an unsolvable
%                   problem are indicated in details instead
%   originalGenes   simply the genes in the input model. Included for
%                   simple presentation of the output
%   details         not all genes will be deleted in all analyses. It is
%                   for example not necessary to delete genes for dead end
%                   reactions. This is a vector with details about
%                   each gene in originalGenes and why or why not it was
%                   deleted
%                   1: Was deleted
%                   2: Proved lethal in sgd (single gene deletion)
%                   3: Only involved in reactions with too many iso-enzymes
%                   4: Involved in dead-end reaction
%
%   NOTE: This function disregards complexes. Any one gene can encode a
%         reaction even if parts of the complex is deleted.
%
%   Usage: [genes, fluxes, originalGenes, details]=findGeneDeletions(model,testType,analysisType,...
%           refModel,oeFactor)
%
%   Rasmus Agren, 2014-01-08
%
originalModel=model;

if nargin<5
    oeFactor=10;
end

%Check that the test type is correct
if ~strcmpi(testType,'sgd') && ~strcmpi(testType,'dgd') && ~strcmpi(testType,'sgo') && ~strcmpi(testType,'dgo')
    EM='Incorrect test type';
    dispEM(EM);
end

%Check that the analysis type is correct
if ~strcmpi(analysisType,'fba') && ~strcmpi(analysisType,'moma')
    EM='Incorrect analysis type';
    dispEM(EM);
end

if (strcmpi(testType,'sgo') || strcmpi(testType,'dgo')) && strcmpi(analysisType,'fba')
    EM='Over expression is only available when using MOMA';
    dispEM(EM);
end

if strcmpi(analysisType,'moma') && nargin<4
    EM='A reference model must be supplied when using MOMA';
    dispEM(EM);
end

originalGenes=model.genes;
details=zeros(numel(model.genes),1);

%First simplify the model to reduce the size
model=simplifyModel(model,true,false,true,true);
model=removeReactions(model,{},true,true); %Removes unused genes
details(~ismember(originalGenes,model.genes))=4;

[~, geneMapping]=ismember(model.genes,originalGenes);

%Get the genes that should be deleted. Not all genes are deleted in all
%optimizations. For example, if all reactions involving a gene have
%iso-enzymes a single deletion of that gene will never have an effect.
%FBA and SGD:   All genes that encode reactions on their own
%FBA and DGD:   All genes that have at most one other gene encoding for at
%               least one of its reactions and doesn't participate
%               exclusively in reactions where a single deletion results in
%               an unsolvable problem
%MOMA and SGD:  All genes that encode reactions on their own
%MOMA and DGD:  All combinations of genes that have at most one other gene
%               encoding for at east one of its reactions and doesn't
%               participate exclusively in reactions where a single
%               deletion results in an unsolvable problem
%MOMA and SGO:  All genes
%MOMA and DGO:  All combinations of genes
genesForSGD=[];
if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
    if strcmpi(testType,'sgd')
        I=sum(model.rxnGeneMat,2)>1; %Reactions with iso-enzymes
    else
        I=sum(model.rxnGeneMat,2)>2; %Reactions with more than one iso-enzyme
    end

    %Find genes involved only in these reactions
    [~, inIsoRxns]=find(model.rxnGeneMat(I,:));
    [~, inNonIsoRxns]=find(model.rxnGeneMat(~I,:));
    genesForSGD=unique(inNonIsoRxns);
    details(geneMapping(setdiff(inIsoRxns,inNonIsoRxns)))=3;
end

%Get the genes that should be deleted in a SGD
if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
   genesToModify=genesForSGD;
end

%Get the genes that should be deleted in SGO
if strcmpi(testType,'sgo')
   genesToModify=1:numel(model.genes);
   genesToModify=genesToModify(:);
end

%Do single deletion/over expression. This is done here since the double
%deletion depends on which single deletions prove lethal
%(to reduce the size of the system)
if strcmpi(testType,'sgd') || strcmpi(testType,'sgo') || strcmpi(testType,'dgd')
    fluxes=zeros(numel(model.rxns),numel(genesToModify));
    solvable=true(numel(genesToModify),1);
    for i=1:numel(genesToModify)
       I=find(model.rxnGeneMat(:,genesToModify(i)));
       if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
           %Constrain all reactions involving the gene to 0
           tempModel=setParam(model,'eq',model.rxns(I),0);
       else
           %To over express a gene, the stoichiometry of the corresponding
           %reactions are changed so that the same flux leads to a higher
           %production
           tempModel=model;
           tempModel.S(:,I)=tempModel.S(:,I).*oeFactor;
       end
       if strcmpi(analysisType,'fba') || strcmpi(testType,'dgd')
            sol=solveLP(tempModel);
       else
            [fluxA, ~, flag]=qMOMA(tempModel,refModel);
            sol.x=fluxA;
            sol.stat=flag;
       end

       %If the optimization terminated successfully
       if sol.stat==1
           fluxes(:,i)=sol.x;
           details(geneMapping(genesToModify(i)))=1;
       else
           solvable(i)=false;
           details(geneMapping(genesToModify(i)))=2;
       end
    end

    fluxes=fluxes(:,solvable);
    genes=geneMapping(genesToModify(solvable));
end

%Now do for DGO. This is rather straight forward since it is always
%solvable and it doesn't matter if there are iso-enzymes
if strcmpi(testType,'dgo')
    genesToModify=nchoosek(1:numel(model.genes),2);
    genes=geneMapping(genesToModify);
    %Since I assume that this is never lethal I set the details already
    details(geneMapping)=1;

    fluxes=sparse(numel(model.rxns),size(genesToModify,1));
    for i=1:size(genesToModify,1)
       I=find(model.rxnGeneMat(:,genesToModify(i,:)));
       %To over express a gene, the stoichiometry of the corresponding
       %reactions are changed so that the same flux leads to a higher
       %production
       tempModel=model;
       tempModel.S(:,I)=tempModel.S(:,I).*oeFactor;
       fluxA=qMOMA(tempModel,refModel);
       fluxes(:,i)=fluxA;
    end
end

%For double deletions FBA or MOMA
if strcmpi(testType,'dgd')
    %This is a little lazy but it's fine. Check which genes that have
    %already been labled as either ony with too many iso-enzymes or
    %non-solveable as a single deletion.
    [~, I]=ismember(originalGenes(details==1),model.genes);
    genesToModify=nchoosek(I,2);
    genes=geneMapping(genesToModify);

    fluxes=sparse(numel(model.rxns),size(genesToModify,1));
    for i=1:size(genesToModify,1)
       I=find(model.rxnGeneMat(:,genesToModify(i,:)));

       %Constrain all reactions involving the gene to 0
       tempModel=setParam(model,'eq',model.rxns(I),0);

       if strcmpi(analysisType,'fba')
            sol=solveLP(tempModel);
       else
            [fluxA, ~, flag]=qMOMA(tempModel,refModel);
            sol.x=fluxA;
            sol.stat=flag;
       end

       if sol.stat==1
           fluxes(:,i)=sol.x;
       end
    end
end

%Map back to the old model
[~, I]=ismember(model.rxns,originalModel.rxns);
temp=fluxes;
fluxes=sparse(numel(originalModel.rxns),size(temp,2));
fluxes(I,:)=temp;
end
