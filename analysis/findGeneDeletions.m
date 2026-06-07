function [genes, fluxes, originalGenes, details, grRatioMuts]=findGeneDeletions(model,testType,analysisType,refModel,oeFactor)
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
%                   (optional, default 'sgd')
%   analysisType    determines whether to use FBA ('fba') or MOMA ('moma')
%                   in the optimization. (optional, default 'fba')
%   refModel        MOMA works by fitting the flux distributions of two
%                   models to be as similar as possible. The most common
%                   application is where you have a reference model where
%                   some of the fluxes are constrained from experimental
%                   data. This model is required when using MOMA
%   oeFactor        a factor by which the fluxes should be increased if a
%                   gene is overexpressed (optional, default 10)
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
%                   1: Was deleted/overexpressed
%                   2: Proved lethal in sgd (single gene deletion)
%                   3: - redundant, no longer used -
%                   4: Involved in dead-end reaction
%   grRatioMuts     growth rate ratio between mutated strain and wild type,
%                   matches the originalGenes(genes) mutants. Note that
%                   this does not directly map to model.genes, as is the case
%                   for COBRA getEssentialGenes. However, this can be 
%                   obtained by afterwards running:
%                       grRatio=zeros(1,numel(model.genes));
%                       grRatio(genes)=grRatioMuts;
%
% Usage: [genes, fluxes, originalGenes, details, grRatioMuts]=findGeneDeletions(model,testType,analysisType,...
%           refModel,oeFactor)

originalModel=model;
if nargin<5
    oeFactor=10;
end
if nargin<2
    testType='sgd';
else
    testType=char(testType);
end

%Check that the test type is correct
if ~strcmpi(testType,'sgd') && ~strcmpi(testType,'dgd') && ~strcmpi(testType,'sgo') && ~strcmpi(testType,'dgo')
    EM='Incorrect test type';
    dispEM(EM);
end

%Check that the analysis type is correct
if nargin<3
    analysisType = 'fba';
else
    analysisType=char(analysisType);
    if ~any(strcmpi(analysisType,{'fba','moma'}))
        EM='Incorrect analysis type';
        dispEM(EM);
    end
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
growthWT=solveLP(model);
growthWT=growthWT.f;

%Do single deletion/over expression. This is done here since the double
%deletion depends on which single deletions prove lethal (to reduce the
%size of the system)
if strcmpi(testType,'sgd') || strcmpi(testType,'sgo') || strcmpi(testType,'dgd')
    fluxes=zeros(numel(model.rxns),numel(model.genes));
    grRatioMuts=zeros(1,numel(model.genes));
    for i=1:numel(model.genes)
        if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
            %Constrain all reactions involving the gene to 0
            tempModel=removeGenes(model,i,false,false,false);
        else
            %To over express a gene, the stoichiometry of the corresponding
            %reactions are changed so that the same flux leads to a higher
            %production
            tempModel=model;
            I=find(model.rxnGeneMat(:,i));
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
            grRatioMuts(i)=sol.f/growthWT;
            details(geneMapping(i))=1;
        else
            fluxes(:,i)=0;
            grRatioMuts(i)=0;
            details(geneMapping(i))=2;
        end
    end
    genes=geneMapping;
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
        grRatioMuts(i)=fluxA(logical(model.c))/growthWT;
    end
end

%For double deletions FBA or MOMA
if strcmpi(testType,'dgd')
    %This is a little lazy but it's fine. Check which genes have already
    %been deleted in 'sgd' analysis.
    [~, I]=ismember(originalGenes(details==1),model.genes);
    genesToModify=nchoosek(I,2);
    genes=geneMapping(genesToModify);
    grRatioMuts=zeros(1,numel(genes));
    fluxes=sparse(numel(model.rxns),size(genesToModify,1));
    for i=1:size(genesToModify,1)
        tempModel=removeGenes(model,genesToModify(i,:),false,false,false);
        
        if strcmpi(analysisType,'fba')
            sol=solveLP(tempModel);
        else
            [fluxA, ~, flag]=qMOMA(tempModel,refModel);
            sol.x=fluxA;
            sol.stat=flag;
        end
        
        if sol.stat==1
            fluxes(:,i)=sol.x;
            grRatioMuts(i)=sol.f/growthWT;
        end
    end
end

%Map back to the old model
[~, I]=ismember(model.rxns,originalModel.rxns);
temp=fluxes;
fluxes=sparse(numel(originalModel.rxns),size(temp,2));
fluxes(I,:)=temp;
end
    