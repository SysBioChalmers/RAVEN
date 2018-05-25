function [rxnScores, geneScores, hpaScores, arrayScores]=scoreModel(model,hpaData,arrayData,tissue,celltype,noGeneScore,multipleGeneScoring,multipleCellScoring,hpaLevelScores)
% scoreRxns
%   Scores the reactions and genes in a model based on expression data
%   from HPA and/or gene arrays
%
%   model               a model structure
%   hpaData             HPA data structure from parseHPA (opt if arrayData is
%                       supplied, default [])
%   arrayData           gene expression data structure (opt if hpaData is
%                       supplied, default [])
%       genes           cell array with the unique gene names
%       tissues         cell array with the tissue names. The list may not be
%                       unique, as there can be multiple cell types per tissue
%       celltypes       cell array with the cell type names for each tissue
%       levels          GENESxTISSUES array with the expression level for
%                       each gene in each tissue/celltype. NaN should be
%                       used when no measurement was performed
%   tissue              tissue to score for. Should exist in either
%                       hpaData.tissues or arrayData.tissues
%   celltype            cell type to score for. Should exist in either
%                       hpaData.celltypes or arrayData.celltypes for this
%                       tissue (opt, default is to use the best values
%                       among all the cell types for the tissue. Use [] if
%                       you want to supply more arguments)
%   noGeneScore         score for reactions without genes (opt, default -2)
%   multipleGeneScoring determines how scores are calculated for reactions
%                       with several genes ('best' or 'average')
%                       (opt, default 'best')
%   multipleCellScoring determines how scores are calculated when several
%                       cell types are used ('best' or 'average')
%                       (opt, default 'best')
%   hpaLevelScores      structure with numerical scores for the expression
%                       level categories from HPA. The structure should have a
%                       "names" and a "scores" field (opt, see code for
%                       default scores)
%
%   rxnScores       scores for each of the reactions in model
%   geneScores      scores for each of the genes in model. Genes which are
%                   not in the dataset(s) have -Inf as scores
%   hpaScores       scores for each of the genes in model if only taking hpaData
%                   into account. Genes which are not in the dataset(s)
%                   have -Inf as scores
%   arrayScores     scores for each of the genes in model if only taking arrayData
%                   into account. Genes which are not in the dataset(s)
%                   have -Inf as scores
%
%   Usage: [rxnScores, geneScores, hpaScores, arrayScores]=scoreModel(model,...
%               hpaData,arrayData,tissue,celltype,noGeneScore,multipleGeneScoring,...
%               multipleCellScoring,hpaLevelScores)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<3
    arrayData=[];
end
if nargin<5
    celltype=[];
end
if nargin<6
    noGeneScore=-2;
end
if nargin<7
    multipleGeneScoring='best';
end
if nargin<8
    multipleCellScoring='best';
end
if nargin<9
    %The first four are for APE, the other ones for staining
    hpaLevelScores.names={'High' 'Medium' 'Low' 'None' 'Strong' 'Moderate' 'Weak' 'Negative'};
    hpaLevelScores.scores=[20 15 10 -8 20 15 10 -8];
end

if isempty(hpaData) && isempty(arrayData)
    EM='Must supply hpaData, arrayData or both';
    dispEM(EM);
end
if ~strcmpi(multipleGeneScoring,'best') && ~strcmpi(multipleGeneScoring,'average')
    EM='Valid options for multipleGeneScoring are "best" or "average"';
    dispEM(EM);
end
if ~strcmpi(multipleCellScoring,'best') && ~strcmpi(multipleCellScoring,'average')
    EM='Valid options for multipleCellScoring are "best" or "average"';
    dispEM(EM);
end


%Throw an error if array data for only one tissue is supplied
if any(arrayData)
    if numel(arrayData.tissues)<2
        EM='arrayData must contain measurements for at least two celltypes/tissues since the score is calculated based on the expression level compared to the overall average';
        dispEM(EM);
    end
end

%This is so that the code can ignore which combination of input data that
%is used
if isempty(arrayData)
    arrayData.genes={};
    arrayData.tissues={};
    arrayData.celltypes={};
    arrayData.levels=[];
end
if isempty(hpaData)
    hpaData.genes={};
    hpaData.tissues={};
    hpaData.celltypes={};
    hpaData.levels={};
    hpaData.types={};
    hpaData.reliabilities={};
    hpaData.gene2Level=[];
    hpaData.gene2Type=[];
    hpaData.gene2Reliability=[];
end

%Check that the tissue exists
if ~ismember(upper(tissue),upper(hpaData.tissues)) && ~ismember(upper(tissue),upper(arrayData.tissues))
    EM='The tissue name does not match';
    dispEM(EM);
end
if any(celltype)
    %Check that both data types has cell type defined if that is to be used
    if ~isfield(hpaData,'celltypes') || ~isfield(arrayData,'celltypes')
        EM='Both hpaData and arrayData must contain cell type information if cell type is to be used';
        dispEM(EM);
    end
    if ~ismember(upper(celltype),upper(hpaData.celltypes)) && ~ismember(upper(celltype),upper(arrayData.celltypes))
        EM='The cell type name does not match';
        dispEM(EM);
    end
end

%Some preprocessing of the structures to speed up a little Remove all
%tissues that are not the correct one
J=~strcmpi(hpaData.tissues,tissue);

%If cell type is supplied, then only keep that cell type
if any(celltype)
    J=J | ~strcmpi(hpaData.celltypes,celltype);
end

hpaData.tissues(J)=[];
if isfield(hpaData,'celltypes')
    hpaData.celltypes(J)=[];
end
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(:,J)=[];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(:,J)=[];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(:,J)=[];
end

%Remove all genes from the structures that are not in model or that aren't
%measured in the tissue
if ~isempty(hpaData.genes) %This should not be necessary, but the summation is a 0x1 matrix and the other is []
    I=~ismember(hpaData.genes,model.genes) | sum(hpaData.gene2Level,2)==0;
else
    I=[];
end
hpaData.genes(I)=[];
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(I,:)=[];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(I,:)=[];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(I,:)=[];
end

I=strcmpi(arrayData.tissues,tissue);
%If cell type is supplied, then only keep that cell type
if any(celltype)
    I=I & strcmpi(arrayData.celltypes,celltype);
end

%Remove all genes from the structures that are not in model or that aren't
%measured in the tissue
J=~ismember(arrayData.genes,model.genes) | myAll(isnan(arrayData.levels(:,I)),2);
arrayData.genes(J)=[];
arrayData.levels(J,:)=[];

%Calculate the scores for the arrayData. These scores are calculated for
%each genes from its fold change between the tissue/celltype(s) in question
%and all other celltypes. This is a lower quality data than protein
%abundance, since a gene that is equally highly expressed in all cell types
%will have a score of 0.0. These scores are therefore only used for genes
%for which there is no HPA data available. The fold changes are transformed
%as min(log(x),10) for x>1 and max(log(x),-5) for x<1 in order to have
%negative scores for lower expressed genes and to scale the scrores to have
%somewhat lower weights than the HPA scores
tempArrayLevels=arrayData.levels;
tempArrayLevels(isnan(tempArrayLevels))=0;
average=sum(tempArrayLevels,2)./sum(~isnan(arrayData.levels),2);
if strcmpi(multipleCellScoring,'best')
    current=max(tempArrayLevels(:,I),[],2);
else
    current=sum(tempArrayLevels(:,I),2)./sum(~isnan(arrayData.levels(:,I)),2);
end
if any(current)
    aScores=5*log(current./average);
else
    aScores=[];
end
aScores(aScores>0)=min(aScores(aScores>0),10);
aScores(aScores<0)=max(aScores(aScores<0),-5);

%Map the HPA levels to scores
[I, J]=ismember(upper(hpaData.levels),upper(hpaLevelScores.names));
if ~all(I)
    EM='There are expression level categories that do not match to hpaLevelScores';
    dispEM(EM);
end
[K, L, M]=find(hpaData.gene2Level);
scores=hpaLevelScores.scores(J);
if strcmpi(multipleCellScoring,'best')
    hScores=max(sparse(K,L,scores(M),numel(hpaData.genes),numel(hpaData.tissues)),[],2);
else
    hScores=mean(sparse(K,L,scores(M),numel(hpaData.genes),numel(hpaData.tissues)),2);
end

%Get the scores for the genes, only use HPA if available
geneScores=inf(numel(model.genes),1)*-1;
hpaScores=geneScores;
arrayScores=geneScores;

[I, J]=ismember(model.genes,hpaData.genes);
hpaScores(I)=hScores(J(I));
geneScores(I)=hScores(J(I));
[I, J]=ismember(model.genes,arrayData.genes);
arrayScores(I)=aScores(J(I));
geneScores(I & myIsInf(geneScores))=aScores(J(I & myIsInf(geneScores)));

%Remove the genes that have no data from the model
I=ismember(model.genes,hpaData.genes) | ismember(model.genes,arrayData.genes);
model.genes(~I)=[];
model.rxnGeneMat(:,~I)=[];

%Map the genes to the HPA/array genes
[hpaExist, hpaMap]=ismember(model.genes,hpaData.genes);
[arrayExist, arrayMap]=ismember(model.genes,arrayData.genes);

%Set the default scores for reactions without genes
rxnScores=ones(numel(model.rxns),1)*noGeneScore;

%Loop through the reactions and calculate the scores
for i=1:numel(model.rxns)
    %Check if it has genes
    I=find(model.rxnGeneMat(i,:));
    if any(I)
        %If any of the genes exist in hpaData, then don't use arrayData
        if any(hpaExist(I))
            %At least one gene was found in HPA
            if strcmpi(multipleGeneScoring,'best')
                rxnScores(i)=max(hScores(hpaMap(I(hpaExist(I)))));
            else
                rxnScores(i)=mean(hScores(hpaMap(I(hpaExist(I)))));
            end
        else
            %Use array data
            if any(arrayExist(I))
                %At least one gene was found in the array data
                if strcmpi(multipleGeneScoring,'best')
                    rxnScores(i)=max(aScores(arrayMap(I(arrayExist(I)))));
                else
                    rxnScores(i)=mean(aScores(arrayMap(I(arrayExist(I)))));
                end
            end
        end
    end
end
end

%This is because isinf and all returns 0x1 for empty set, which gives a
%concatenation error. Do like this instead of having many if statements
function y=myIsInf(x)
y=isinf(x);
if isempty(y)
    y=[];
end
end
function y=myAll(x,dim)
y=all(x,dim);
if isempty(y)
    y=[];
end
end
