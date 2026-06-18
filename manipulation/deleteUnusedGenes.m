function reducedModel=deleteUnusedGenes(model,varargin)
% deleteUnusedGenes  Delete all genes not associated to any reaction.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% verbose : double
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

if verbose && isfield(model,'rxnGeneMat') && isfield(model,'genes')
    [~, b]=find(model.rxnGeneMat);
    toKeep=false(numel(model.genes),1);
    toKeep(b)=true;
    switch verbose
        case 1
            fprintf('Number of unused genes removed from the model: %d\n', sum(~toKeep))
        case 2
            disp('The following genes were removed from the model:')
            disp(model.genes(~toKeep))
    end
end

reducedModel=removeReactions(model,[],'removeUnusedGenes',true);
end
