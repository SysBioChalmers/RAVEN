function reducedModel=deleteUnusedGenes(model,verbose)
% deleteUnusedGenes
%   Deletes all genes that are not associated to any reaction
%
%   model           a model structure
%   verbose         0 for silent; 1 for printing number of deleted genes;
%                   2 for printing the list of deleted genes (optional, default 1)
%
%   reducedModel    an updated model structure
%
% Usage: reducedModel=deleteUnusedGenes(model)

if nargin<2
    verbose=1;
end
reducedModel=model;

%Find all genes that are not used
[~, b]=find(reducedModel.rxnGeneMat);
toKeep=false(numel(reducedModel.genes),1);
toKeep(b)=true;

switch verbose
    case 1
        disp('Number of unused genes removed from the model:')
        disp(numel(toKeep(~toKeep)))
    case 2
        disp('The following genes were removed from the model:')
        disp(reducedModel.genes(~toKeep))
    case 0
end
        
reducedModel.genes=reducedModel.genes(toKeep);
reducedModel.rxnGeneMat=reducedModel.rxnGeneMat(:,toKeep);

if isfield(reducedModel,'geneShortNames')
    reducedModel.geneShortNames=reducedModel.geneShortNames(toKeep);
end

if isfield(reducedModel,'proteins')
    reducedModel.proteins=reducedModel.proteins(toKeep);
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
