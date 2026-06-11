function reducedModel=deleteUnusedGenes(model,varargin)
% deleteUnusedGenes  Delete all genes not associated to any reaction.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% verbose : double, optional
%     0 for silent; 1 for printing the number of deleted genes; 2 for
%     printing the list of deleted genes (default 1).
%
% Returns
% -------
% reducedModel : struct
%     an updated model structure.
%
% Examples
% --------
%     reducedModel=deleteUnusedGenes(model);

p=parseRAVENargs(varargin, {'verbose',1});
verbose=p.verbose;
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
