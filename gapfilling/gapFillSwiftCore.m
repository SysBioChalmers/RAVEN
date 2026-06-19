function activeRxns = gapFillSwiftCore(model, coreRxns, epsilon)
% gapFillSwiftCore  SWIFTCORE LP subroutine for swiftGapFill.
%
% Single-LP alternative to gapFillFastCore. Maximises the sum of non-core
% fluxes while enforcing that every core reaction carries at least epsilon
% flux. Any reaction with non-zero flux in the solution is "consistent"
% with the core.
%
% Because the maximisation objective makes the LP degenerate (many optima
% exist that differ only in which non-core reactions carry flux), the
% returned active set can vary between LP solvers and between runs with
% slightly different numerical conditions. This is the expected behaviour
% of SWIFTCORE (Tefagh & Boyd 2020).
%
% Based on SWIFTCORE (Tefagh & Boyd 2020, BMC Bioinformatics 21:23).
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
% coreRxns : cell or logical or double
%     core reaction set: cell array of reaction IDs, a logical vector, or a
%     vector of reaction indices into model.rxns. All core reactions are
%     forced to carry flux >= epsilon in the LP solution.
% epsilon : double
%     minimum flux threshold for a reaction to be considered active.
%
% Returns
% -------
% activeRxns : logical
%     logical vector (length numel(model.rxns)) where true indicates the
%     reaction carries flux in the solution. Empty on infeasibility.
%
% See also
% --------
% gapFillFastCore, gapFillFastLP

if nargin < 3 || isempty(epsilon)
    epsilon = 1e-4;
end

% ---- Resolve core reaction indices ----
% coreOrigIdx: numeric indices into model.rxns
if iscell(coreRxns)
    coreOrigIdx = find(getIndexes(model, coreRxns, 'rxns', true));
elseif islogical(coreRxns)
    coreOrigIdx = find(coreRxns);
else
    coreOrigIdx = coreRxns(:);
end

% ---- Convert to irreversible form ----
% SWIFTCORE requires an irreversible model so all fluxes are non-negative
% and the maximisation objective is well-defined.
[irrevModel, ~, ~, irrev2rev] = convertToIrrev(model);
nIrrev = numel(irrevModel.rxns);
nMets  = numel(irrevModel.mets);

% Map core to irreversible equivalents.
coreIdxIrrev = ismember(irrev2rev, coreOrigIdx);

% ---- Formulate LP ----
% Maximise sum of non-core fluxes, forcing core reactions >= epsilon.
%
%   max   sum_{j not in C} v_j
%   s.t.  S_irrev * v = 0
%         v_j >= epsilon   for j in C
%         0 <= v <= ub
%
lb = irrevModel.lb;
lb(coreIdxIrrev) = max(lb(coreIdxIrrev), epsilon);

c = -ones(nIrrev, 1);    % negate because optimizeProb minimises by default
c(coreIdxIrrev) = 0;     % core reactions not in objective

prob.A      = irrevModel.S;
prob.a      = prob.A;
prob.b      = zeros(nMets, 1);
prob.csense = repmat('E', 1, nMets);
prob.c      = c;
prob.osense = 1;           % minimise (with negated c → maximises the original)
prob.lb     = lb;
prob.ub     = irrevModel.ub;
prob.vartype = repmat('C', 1, nIrrev);

sol = optimizeProb(prob, [], false);

if ~checkSolution(sol)
    activeRxns = false(numel(model.rxns), 1);
    return;
end

% ---- Map solution back to original model ----
irrevActive = sol.full(1:nIrrev) >= epsilon / 2;
activeRxns  = false(numel(model.rxns), 1);
for k = 1:nIrrev
    if irrevActive(k)
        activeRxns(irrev2rev(k)) = true;
    end
end
end
