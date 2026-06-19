function [solutions, info] = sampleFluxSpace(model, varargin)
% sampleFluxSpace  Near-uniform MCMC sampling of the flux solution space.
%
% Draws flux vectors approximately uniformly from the flux polytope
% {v : S*v = 0, lb <= v <= ub} by one of two Markov-chain Monte Carlo methods:
%
%   'chrr' (default) — Coordinate Hit-and-Run with Rounding (Haraldsdottir et al.
%       2017). The polytope is rounded with its maximum-volume inscribed
%       ellipsoid before sampling, so mixing is independent of how elongated the
%       feasible set is. Recommended for enzyme-constrained (ecModel +
%       proteomics) and flux-measured models, whose feasible set is a thin,
%       ill-conditioned slab that defeats unrounded chains.
%   'achr' — Artificially Centered Hit-and-Run (Kaufman & Smith 1998). Lighter on
%       well-conditioned models (no rounding), but slower-mixing on elongated
%       polytopes.
%
% This complements randomSampling, which returns diverse *vertices* of the
% polytope (random-objective method); sampleFluxSpace returns the (near-)uniform
% *interior* distribution. Set any constraints to condition on (e.g. a biomass
% lower bound, measured exchange fluxes) on the model before calling.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
%
% Name-Value Arguments
% --------------------
% method : char
%     'chrr' (default) or 'achr'.
% nSamples : double
%     number of flux vectors to return (default 1000).
% thinning : double
%     Markov-chain steps between recorded samples (default 100).
% nBurnin : double
%     CHRR only: burn-in steps before the first recorded sample (default 1000).
% seed : double
%     random seed for reproducibility (default: current RNG state).
% fixedWidthTol : double
%     CHRR only: flux-range threshold below which a reaction is folded into the
%     equality system as fixed (default 1e-7).
%
% Returns
% -------
% solutions : double
%     nRxns-by-nSamples matrix of sampled flux distributions.
% info : struct
%     method-specific diagnostics. For 'chrr': nDimensions, mveConverged,
%     fixedRxns. For 'achr': nWarmup (number of warmup vertices).
%
% Examples
% --------
%     sols = sampleFluxSpace(model, 'method', 'chrr', 'nSamples', 500, 'seed', 1);
%     sols = sampleFluxSpace(model, 'method', 'achr', 'nSamples', 500);
%
% See also
% --------
% sampleCHRR, sampleACHR, sampleMaxVolEllipse, randomSampling

p = parseRAVENargs(varargin, { ...
    'method',        'chrr'; ...
    'nSamples',      1000; ...
    'thinning',      100; ...
    'nBurnin',       1000; ...
    'seed',          []; ...
    'fixedWidthTol', 1e-7});

if p.nSamples <= 0
    error('RAVEN:badInput', 'sampleFluxSpace: nSamples must be positive.');
end

switch lower(p.method)
    case 'chrr'
        [solutions, info] = sampleCHRR(model, p.nSamples, p.thinning, ...
            p.nBurnin, p.seed, p.fixedWidthTol);
        info.method = 'chrr';
    case 'achr'
        warmupPts = sampleWarmupPoints(model);
        solutions = sampleACHR(model, p.nSamples, p.thinning, warmupPts, p.seed);
        info = struct('method', 'achr', 'nWarmup', size(warmupPts, 1));
    otherwise
        error('RAVEN:badInput', ...
            'sampleFluxSpace: unknown method ''%s''; expected ''chrr'' or ''achr''.', ...
            p.method);
end
end
