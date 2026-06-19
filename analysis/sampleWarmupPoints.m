function warmup = sampleWarmupPoints(model)
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
% Returns
% -------
% warmup : double
%     nPoints-by-nRxns matrix; each row is a vertex flux distribution. Duplicate
%     rows are removed.
%
% See also
% --------
% sampleACHR, randomSampling

nRxns = numel(model.rxns);
pts = zeros(2*nRxns, nRxns);
np = 0;

tmp = model;
for d = [1, -1]                        % maximise, then minimise each reaction
    for i = 1:nRxns
        if model.ub(i) - model.lb(i) < 1e-9
            continue;                  % fixed reaction — no warmup direction
        end
        tmp.c = zeros(nRxns, 1);
        tmp.c(i) = d;                  % solveLP maximises c'*v
        sol = solveLP(tmp);
        if isempty(sol.x) || sol.stat <= 0
            continue;
        end
        np = np + 1;
        pts(np, :) = sol.x';
    end
end

warmup = pts(1:np, :);
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
