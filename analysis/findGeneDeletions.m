function [genes, fluxes, originalGenes, details, grRatioMuts]=findGeneDeletions(model,varargin)
% findGeneDeletions  Delete genes and track the resulting fluxes.
%
% Deletes genes, optimizes the model by FBA, and keeps track of the
% resulting fluxes. This is used for identifying gene deletion targets.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% testType : char
%     single/double gene deletion (default "sgd"):
%
%     - "sgd" : single gene deletion
%     - "dgd" : double gene deletion
%
% Returns
% -------
% genes : double
%     a matrix with the genes that were deleted in each optimization (the
%     gene indexes in originalGenes). Each row corresponds to a column in
%     fluxes.
% fluxes : double
%     a matrix with the resulting fluxes. Double deletions that result in
%     an unsolvable problem have all zero flux. Single deletions that
%     result in an unsolvable problem are indicated in details instead.
% originalGenes : cell
%     simply the genes in the input model. Included for simple
%     presentation of the output.
% details : double
%     not all genes will be deleted in all analyses. It is for example not
%     necessary to delete genes for dead end reactions. This is a vector
%     with details about each gene in originalGenes and why or why not it
%     was deleted:
%
%     - 1 : was deleted
%     - 2 : proved lethal in sgd (single gene deletion)
%     - 3 : redundant, no longer used
%     - 4 : involved in dead-end reaction
% grRatioMuts : double
%     growth rate ratio between mutated strain and wild type, matching the
%     originalGenes(genes) mutants. Note that this does not directly map
%     to model.genes, as is the case for COBRA getEssentialGenes. However,
%     this can be obtained by afterwards running:
%
%         grRatio=zeros(1,numel(model.genes));
%         grRatio(genes)=grRatioMuts;
%
% Examples
% --------
%     [genes, fluxes, originalGenes, details, grRatioMuts]=...
%         findGeneDeletions(model,testType);

originalModel=model;
p=parseRAVENargs(varargin, {'testType',[]});
testType=p.testType;
if isempty(testType)
    testType='sgd';
else
    testType=char(testType);
end

%Check that the test type is correct
if ~strcmpi(testType,'sgd') && ~strcmpi(testType,'dgd')
    EM='Incorrect test type';
    error('RAVEN:badInput', '%s', EM);
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

%Do single deletion. This is done here since the double deletion depends
%on which single deletions prove lethal (to reduce the size of the system)
if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
    fluxes=zeros(numel(model.rxns),numel(model.genes));
    grRatioMuts=zeros(1,numel(model.genes));
    for i=1:numel(model.genes)
        %Constrain all reactions involving the gene to 0
        tempModel=removeGenes(model,i,false,false,false);
        sol=solveLP(tempModel);

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

%For double deletions
if strcmpi(testType,'dgd')
    %This is a little lazy but it is fine. Check which genes have already
    %been deleted in "sgd" analysis.
    [~, I]=ismember(originalGenes(details==1),model.genes);
    genesToModify=nchoosek(I,2);
    genes=geneMapping(genesToModify);
    grRatioMuts=zeros(1,numel(genes));
    fluxes=sparse(numel(model.rxns),size(genesToModify,1));
    for i=1:size(genesToModify,1)
        tempModel=removeGenes(model,genesToModify(i,:),false,false,false);
        sol=solveLP(tempModel);

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
