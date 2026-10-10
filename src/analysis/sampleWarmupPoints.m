function warmup = sampleWarmupPoints(model,varargin)
% sampleWarmupPoints  FVA warmup vertices for ACHR sampling.
%
% Generates warmup points for sampleACHR by maximising and minimising each
% reaction in turn and storing the full flux distribution at each optimum.
% These vertices span the flux polytope and seed the hit-and-run directions.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
%
% Name-Value Arguments
% --------------------
% runParallel : logical
%     speed up calculations by parallel processing, as in getAllowedBounds,
%     which solves the same per-reaction min/max LPs (default true).
%
% Returns
% -------
% warmup : double
%     nPoints-by-nRxns matrix; each row is a vertex flux distribution. Duplicate
%     rows are removed.
%
% See also
% --------
% sampleACHR, randomSampling, getAllowedBounds

p=parseRAVENargs(varargin, {'runParallel',true});
runParallel=p.runParallel;

nRxns = numel(model.rxns);
nW = parallelWorkersRAVEN(runParallel);
fixedRxn = (model.ub - model.lb) < 1e-9; % fixed reaction — no warmup direction

%Maximise and minimise each reaction in turn. Each direction is its own
%fixed-size (rxns x rxns) buffer (NaN where a reaction is fixed or
%infeasible) so both parfor loops write to statically-indexed rows,
%mirroring getAllowedBounds' own min/max-per-reaction parfor.
ptsMax = NaN(nRxns, nRxns);
ptsMin = NaN(nRxns, nRxns);

PB = progressReport(2*nRxns,'Running sampleWarmupPoints');
parfor (i = 1:nRxns, nW)
    if ~fixedRxn(i)
        tmp = model;
        tmp.c = zeros(nRxns, 1);
        tmp.c(i) = 1;                  % solveLP maximises c'*v
        sol = solveLP(tmp);
        if ~isempty(sol.x) && sol.stat > 0
            ptsMax(i,:) = sol.x';
        end
    end
    count(PB)
end
parfor (i = 1:nRxns, nW)
    if ~fixedRxn(i)
        tmp = model;
        tmp.c = zeros(nRxns, 1);
        tmp.c(i) = -1;
        sol = solveLP(tmp);
        if ~isempty(sol.x) && sol.stat > 0
            ptsMin(i,:) = sol.x';
        end
    end
    count(PB)
end

pts = [ptsMax; ptsMin];
warmup = pts(~any(isnan(pts),2), :);
if isempty(warmup)
    error('RAVEN:sampling', ...
        'sampleWarmupPoints: could not generate any warmup points; check model feasibility.');
end
warmup = unique(warmup, 'rows');
if size(warmup, 1) < 2
    error('RAVEN:sampling', ...
        'sampleWarmupPoints: flux cone collapses to a single point; cannot sample.');
end
end
