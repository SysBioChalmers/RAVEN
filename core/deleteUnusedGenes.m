function reducedModel=deleteUnusedGenes(model)
% deleteUnusedGenes
%   Deletes all genes that are not associated to any reaction
%
%   model           a model structure
%
%   reducedModel    an updated model structure
%
%   Usage: reducedModel=deleteUnusedGenes(model)
%
%   Simonas Marcisauskas, 2017-09-12
%

reducedModel=model;

%Find all genes that are not used
[~, b]=find(reducedModel.rxnGeneMat);
toKeep=false(numel(reducedModel.genes),1);
toKeep(b)=true;

reducedModel.genes=reducedModel.genes(toKeep);
reducedModel.rxnGeneMat=reducedModel.rxnGeneMat(:,toKeep);

disp('Number of unused genes removed from the model:')
disp(numel(toKeep(~toKeep)))

if isfield(reducedModel,'geneShortNames')
    reducedModel.geneShortNames=reducedModel.geneShortNames(toKeep);
end

if isfield(reducedModel,'geneMiriams')
    reducedModel.geneMiriams=reducedModel.geneMiriams(toKeep);
end

if isfield(reducedModel,'geneFrom')
    reducedModel.geneFrom=reducedModel.geneFrom(toKeep);
end

if isfield(reducedModel,'geneComps')
    reducedModel.geneComps=reducedModel.geneComps(toKeep);
end
end
