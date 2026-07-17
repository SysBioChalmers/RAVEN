function [solutions, info] = sampleCHRR(model, nSamples, thinning, nBurnin, seed, fixedWidthTol)
% sampleCHRR  Coordinate Hit-and-Run with Rounding flux sampling.
%
% Draws nSamples flux vectors approximately uniformly from the flux polytope
% {v : S*v = 0, lb <= v <= ub}, using the CHRR algorithm (Haraldsdottir et al.
% 2017, Bioinformatics 33:1741):
%
%   1. Reduce the polytope to a full-dimensional body via the nullspace of S
%      (v = v0 + N*x); reactions with ~zero flux range are folded into the
%      equality system so the reduced polytope is full-dimensional.
%   2. Round it with the maximum-volume inscribed ellipsoid (sampleMaxVolEllipse).
%      Rounding makes mixing independent of how elongated the polytope is, which
%      is the regime of enzyme-constrained (ecModel + proteomics) and
%      flux-measured models.
%   3. Walk the rounded polytope with coordinate hit-and-run and map back.
%
% Set any constraints you want to condition on (e.g. a biomass lower bound or
% measured exchange fluxes) on the model before calling.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
% nSamples : double
%     number of flux vectors to return (default 1000).
% thinning : double
%     Markov-chain steps between recorded samples (default 100).
% nBurnin : double
%     burn-in steps discarded before the first recorded sample (default 1000).
% seed : double
%     random seed for reproducibility (default: current RNG state).
% fixedWidthTol : double
%     a reaction whose allowed flux range is narrower than this is treated as
%     stoichiometrically fixed and folded into the equality system, keeping the
%     reduced polytope full-dimensional (default 1e-7).
%
% Returns
% -------
% solutions : double
%     nRxns-by-nSamples matrix of sampled flux distributions.
% info : struct with fields:
%   .nDimensions  — dimension of the sampled flux polytope
%   .mveConverged — whether MVE rounding reached tolerance
%   .fixedRxns    — IDs of reactions folded in as fixed
%
% See also
% --------
% randomSampling, sampleACHR, sampleMaxVolEllipse

if nargin < 2 || isempty(nSamples);      nSamples = 1000;     end
if nargin < 3 || isempty(thinning);      thinning = 100;      end
if nargin < 4 || isempty(nBurnin);       nBurnin = 1000;      end
if nargin < 5;                           seed = [];           end
if nargin < 6 || isempty(fixedWidthTol); fixedWidthTol = 1e-7; end
if ~isempty(seed)
    rng(seed);
end

nRxns = numel(model.rxns);
S  = full(model.S);
lb = model.lb(:);
ub = model.ub(:);

% ---- Fold implicitly-fixed reactions into equalities ----
% A reaction with ~zero flux range is determined by stoichiometry; keeping it as
% a box constraint would make the reduced polytope lower-dimensional and the MVE
% solve singular. Fix it at its (unique) value instead.
[minF, maxF] = getAllowedBounds(model);
minF = minF(:); maxF = maxF(:);
width = maxF - minF;
fixed = width < fixedWidthTol;
fixedVals = (maxF + minF) / 2;
fixedIdx = find(fixed);

Aeq = S;
beq = zeros(size(S, 1), 1);
if ~isempty(fixedIdx)
    blk = zeros(numel(fixedIdx), nRxns);
    for k = 1:numel(fixedIdx)
        blk(k, fixedIdx(k)) = 1;
    end
    Aeq = [Aeq; blk];
    beq = [beq; fixedVals(fixedIdx)];
end

% Particular solution v0 and nullspace basis N: v = v0 + N*x.
% v0 is ANY particular solution; the sampled flux distribution is invariant to
% its choice and to the nullspace basis, so MATLAB's '\' (a sparse basic
% solution) and the Python reference's lstsq (minimum-norm solution) sample the
% same polytope despite returning different v0.
v0 = Aeq \ beq;
N  = null(Aeq);
d  = size(N, 2);

info = struct();
info.fixedRxns = model.rxns(fixedIdx);

if d == 0
    % Flux space is a single point — nothing to sample.
    solutions = repmat(v0, 1, nSamples);
    info.nDimensions = 0;
    info.mveConverged = true;
    return;
end

% ---- Build the full-dimensional polytope in x (free reactions only) ----
free = ~fixed;
Nf = N(free, :);
A_full = [Nf; -Nf];
b_full = [ub(free) - v0(free); v0(free) - lb(free)];
% Drop rows whose direction vanished in the nullspace (constraint is constant).
keep = sqrt(sum(A_full.^2, 2)) > 1e-10;
A_full = A_full(keep, :);
b_full = b_full(keep);

% ---- Round with the MVE ----
x0 = sampleChebyshevCenter(A_full, b_full);
[center, E, converged] = sampleMaxVolEllipse(A_full, b_full, x0);
info.nDimensions = d;
info.mveConverged = converged;
if ~converged
    warning('RAVEN:warning', '%s', ['The maximum-volume ellipsoid rounding ' ...
        'did not converge; the samples may be poorly mixed. Inspect the ' ...
        'mveConverged field of the second output.']);
end

% Rounded polytope {y : A_r*y <= b_r}, which contains the unit ball.
A_r = A_full * E;
b_r = b_full - A_full*center;

% ---- Coordinate hit-and-run on the rounded polytope ----
y = zeros(d, 1);
s = b_r - A_r*y;                       % maintained slack, > 0 at start
solutions = zeros(nRxns, nSamples);
rec = 0;
total = nBurnin + thinning*nSamples;
colTol = 1e-12;

for stepno = 1:total
    j = randi(d);
    col = A_r(:, j);
    pos = col > colTol;
    neg = col < -colTol;
    ratios = s ./ col;
    if any(pos); tHi = min(ratios(pos)); else; tHi = inf;  end
    if any(neg); tLo = max(ratios(neg)); else; tLo = -inf; end
    if ~isfinite(tHi) || ~isfinite(tLo) || tHi <= tLo
        continue;                      % unbounded/degenerate coordinate; skip
    end
    t = tLo + rand*(tHi - tLo);
    y(j) = y(j) + t;
    s = s - col*t;

    if stepno > nBurnin && mod(stepno - nBurnin, thinning) == 0
        rec = rec + 1;
        x = center + E*y;
        solutions(:, rec) = v0 + N*x;
    end
end

if rec < nSamples
    if rec >= 1
        solutions(:, rec+1:end) = repmat(solutions(:, rec), 1, nSamples - rec);
    else
        solutions = repmat(v0, 1, nSamples);
    end
end
end
