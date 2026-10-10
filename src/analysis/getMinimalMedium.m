function [medium, mediumIdx] = getMinimalMedium(model, varargin)
% getMinimalMedium  Find the minimal set of uptake reactions for growth.
%
% Uses MILP to find the smallest set of exchange reactions that can carry
% uptake flux (lb < 0) while still allowing the model to achieve a target
% growth rate. Unlike iterative heuristics, the MILP formulation finds the
% globally minimal medium.
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure.
%
% Name-Value Arguments
% --------------------
% minGrowth : double
%     minimum objective value (growth rate) the medium must support.
%     Defaults to 10 % of the unconstrained optimum computed automatically.
%     Set explicitly to avoid an extra FBA call (e.g. 'minGrowth', 0.1).
% verbose : logical
%     print the resulting minimal medium (default true).
%
% Returns
% -------
% medium : cell
%     cell array of reaction IDs that form the minimal medium. Empty if
%     no feasible medium can support the requested growth.
% mediumIdx : double
%     column vector of reaction indices corresponding to medium.
%
% Examples
% --------
%     [med, idx] = getMinimalMedium(model)
%     [med, idx] = getMinimalMedium(model, 'minGrowth', 0.05)
%     [med, idx] = getMinimalMedium(model, 'verbose', false)
%
% Notes
% -----
% Only exchange reactions with lb < 0 are considered as candidate uptake
% reactions. Exchange reactions that are always secretion-only (lb >= 0) or
% that are completely blocked are not included.
%
% The MILP formulation:
%   minimise  sum(y_i)
%   s.t.      S v = 0
%             c' v >= minGrowth
%             v_i >= lb_i * y_i   for each candidate uptake reaction i
%             lb_j <= v_j <= ub_j for all reactions j
%             y_i in {0, 1}
%
% A y_i = 0 forces reaction i to carry non-negative flux (no net uptake);
% y_i = 1 allows it to carry flux as low as lb_i.

p = parseRAVENargs(varargin, {'minGrowth',[]; 'verbose',true});
minGrowth = p.minGrowth;
verbose   = p.verbose;

medium    = {};
mediumIdx = zeros(0,1);

% ---- Resolve minGrowth ----
if isempty(minGrowth)
    sol0 = solveLP(model);
    if sol0.stat ~= 1
        error('RAVEN:infeasible', ...
            'getMinimalMedium: model is infeasible under current bounds.');
    end
    minGrowth = 0.1 * sol0.f;
end

if minGrowth <= 0
    warning('getMinimalMedium:zeroGrowth', ...
        'minGrowth <= 0; the MILP will accept a zero-growth solution.');
end

% ---- Identify candidate uptake reactions ----
[~, exchIdx] = getExchangeRxns(model);
% Only consider exchange reactions with lb < 0 (able to take up)
lb_cap = max(model.lb(exchIdx), -1000);  % cap at -1000 for numerical safety
isUptake = lb_cap < 0;
uptakeIdx = exchIdx(isUptake);           % indices into model.rxns
uptakeLb  = lb_cap(isUptake);           % effective lower bounds (negative)
n_uptake  = numel(uptakeIdx);

if n_uptake == 0
    if verbose
        fprintf('getMinimalMedium: no exchange reactions with lb < 0 found.\n');
    end
    return
end

% ---- Build MILP ----
n_mets = numel(model.mets);
n_rxns = numel(model.rxns);
n_vars = n_rxns + n_uptake;   % [flux vars | binary uptake indicators]

% Objective: minimise sum(y)
obj = [zeros(n_rxns,1); ones(n_uptake,1)];

% 1. Steady-state: S*v = 0   (n_mets rows)
A_ss = [model.S, sparse(n_mets, n_uptake)];
b_ss = zeros(n_mets, 1);

% 2. Growth lower bound: c'*v >= minGrowth   (1 row)
A_grw = [model.c', sparse(1, n_uptake)];
b_grw = minGrowth;

% 3. Uptake coupling: v_i - lb_i*y_k >= 0   (n_uptake rows)
%    coefficient of v_i (in flux block) = 1
%    coefficient of y_k (in binary block, k-th uptake var) = -lb_i (> 0)
rows  = (1:n_uptake)';
colsV = uptakeIdx;
colsY = n_rxns + (1:n_uptake)';
A_upt = sparse([rows; rows], [colsV; colsY], [ones(n_uptake,1); -uptakeLb], n_uptake, n_vars);
b_upt = zeros(n_uptake, 1);

% Assemble — prob.a is the "original" A (without solver slack columns) and
% is used by optimizeProb to know how many variables to return in res.full.
prob.A       = [A_ss; A_grw; A_upt];
prob.a       = prob.A;
prob.b       = [b_ss; b_grw; b_upt];
prob.csense  = [repmat('E',1,n_mets), 'G', repmat('G',1,n_uptake)];
prob.c       = obj;
prob.osense  = 1;   % minimise
lb_forProb              = model.lb;
lb_forProb(uptakeIdx)   = uptakeLb;   % cap at -1000 for numerical safety
prob.lb      = [lb_forProb; zeros(n_uptake,1)];
prob.ub      = [model.ub; ones(n_uptake,1)];
prob.vartype = [repmat('C', 1, n_rxns), repmat('I', 1, n_uptake)];

% ---- Solve ----
sol = optimizeProb(prob, [], false);
isFeasible = checkSolution(sol);

if ~isFeasible
    if verbose
        fprintf('getMinimalMedium: MILP returned no feasible solution.\n');
        fprintf('  The model may not be able to grow at minGrowth=%.4g with any medium.\n', minGrowth);
    end
    return
end

% ---- Extract result ----
y_sol     = sol.full(n_rxns + 1 : end);
inMedium  = y_sol > 0.5;
mediumIdx = uptakeIdx(inMedium);
medium    = model.rxns(mediumIdx);

if verbose
    printMedium_(model, medium, mediumIdx, minGrowth, n_uptake);
end
end

%--------------------------------------------------------------------------
function printMedium_(model, medium, mediumIdx, minGrowth, nCandidates)
    W = 62;
    fprintf('\n%s\n', repmat('=', 1, W));
    fprintf('  Minimal medium  (%d of %d candidate uptake reactions)\n', ...
        numel(medium), nCandidates);
    fprintf('  Target growth: %.4g\n', minGrowth);
    fprintf('%s\n', repmat('-', 1, W));
    if isempty(medium)
        fprintf('  (none — model can grow without any uptake)\n');
    else
        for i = 1:numel(medium)
            ri  = mediumIdx(i);
            nm  = '';
            if isfield(model,'rxnNames') && ri <= numel(model.rxnNames)
                nm = char(model.rxnNames{ri});
            end
            if numel(nm) > 35, nm = [nm(1:32) '...']; end
            fprintf('  %-20s  %s\n', medium{i}, nm);
        end
    end
    fprintf('%s\n\n', repmat('=', 1, W));
end
