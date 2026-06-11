function out = runSimpleOptKnock(model, targetRxn, biomassRxn, varargin)
% runSimpleOptKnock  Simple OptKnock for growth-coupled production.
%
% Simple OptKnock algorithm that checks all gene or reaction deletions for
% growth-coupled metabolite production, by testing all possible
% combinations. This is not defined as MILP, and is therefore slow (but
% simple).
%
% Parameters
% ----------
% model : struct
%     a model structure.
% targetRxn : char
%     identifier of target reaction.
% biomassRxn : char
%     identifier of biomass reaction.
%
% Name-Value Arguments
% --------------------
% deletions : cell
%     cell array with gene or reaction identifiers that should be
%     considered for knockout (default model.rxns).
% genesOrRxns : char
%     string indicating whether deletions parameter is given with 'genes'
%     or 'rxns' identifiers (default 'rxns').
% maxNumKO : double
%     maximum number of simultaneous knockouts (default 1).
% minGrowth : double
%     minimum growth rate (default 0.05).
%
% Returns
% -------
% out : struct
%     structure with deletion strategies that result in growth-coupled
%     production, with fields:
%
%     - KO : cell array with gene(s) or reaction(s) to be deleted.
%     - growthRate : vector with growth rates after deletion.
%     - prodRate : vector with production rates after deletion.
%
% Examples
% --------
%     out = runSimpleOptKnock(model, targetRxn, biomassRxn, deletions, ...
%         genesOrRxns, maxNumKO, minGrowth);

p=parseRAVENargs(varargin, {'deletions',[]; 'genesOrRxns','rxns'; 'maxNumKO',1; 'minGrowth',0.05});
deletions=p.deletions;
genesOrRxns=p.genesOrRxns;
maxNumKO=p.maxNumKO;
minGrowth=p.minGrowth;
if isempty(deletions)
    params.deletions = model.rxns;
else
    params.deletions = deletions;
end
params.genesOrRxns = genesOrRxns;
params.maxNumKO = maxNumKO;
params.minGrowth = minGrowth;

% Number of deletions
out.KO         = cell(0,params.maxNumKO); % The KO genes/rxns
out.growthRate = zeros(0);
out.prodRate   = zeros(0);
out.score      = zeros(0);

params.biomassIdx  = getIndexes(model,biomassRxn,'rxns');
params.targetIdx   = getIndexes(model,targetRxn,'rxns');

model = setParam(model,'obj',params.biomassIdx,1);
[solWT, hsSol] = solveLP(model);
WT.minScore = solWT.x(params.targetIdx)*solWT.f;

fprintf('Running simple OptKnock analysis...   0%% complete');
KO=zeros(1,params.maxNumKO);
[~,~,out,~] = knockoutIteration(model,params,WT,out,params.maxNumKO,KO,[],hsSol);

if size(out.KO,2)>1
    singleKO = cellfun(@isempty,out.KO(:,1));
    if any(singleKO)
        singleKO = out.KO{singleKO,2};
        singleKO = strcmp(out.KO(:,1),singleKO);
        out.KO(singleKO,:) = [];
        out.growthRate(singleKO,:) = [];
        out.prodRate(singleKO,:) = [];
        out.score(singleKO,:) = [];
    end
end


fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');
end

function [model,iteration,out,KO,hsSol] = knockoutIteration(model,params,WT,out,iteration,KO,minScore,hsSol)
if nargin<7 || isempty(minScore)
    minScore = WT.minScore;
end
iteration = iteration - 1;
for i = 1:numel(params.deletions)
    if iteration+1==params.maxNumKO
        progress=pad(num2str(floor(i/numel(params.deletions)*100)),3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
    KO(iteration+1)=i;
    switch params.genesOrRxns
        case 'rxns'
            modelKO = setParam(model,'eq',params.deletions{i},0);
        case 'genes'
            modelKO = removeGenes(model,params.deletions{i},false,false,false);
    end
    solKO = solveLP(modelKO,0,[],hsSol);
    if ~isempty(solKO.f)
        growthRate = solKO.f;
        prodRate   = solKO.x(params.targetIdx);
        prodRate(prodRate<1e-10)=0; % Filter out results from solver tolerance
        if growthRate > params.minGrowth
            if growthRate*prodRate > minScore*1.01
                out.KO(end+1,find(KO)) = transpose(params.deletions(KO(find(KO))));
                out.growthRate(end+1,1) = growthRate;
                out.prodRate(end+1,1)   = prodRate;
                out.score(end+1,1)      = growthRate*prodRate;
            end
            if iteration>0
                [~,~,out] = knockoutIteration(modelKO,params,WT,out,iteration,KO,growthRate*prodRate,hsSol);
            end
        end
    end
end
end
