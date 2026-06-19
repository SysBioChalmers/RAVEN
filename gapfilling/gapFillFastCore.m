function activeRxns = gapFillFastCore(model, coreRxns, epsilon)
% gapFillFastCore  L1-norm LP subroutine for fastGapFill.
%
% Finds the minimum-L1-flux extension of the model that activates all
% reactions in the core set. Any reaction with non-zero flux in the
% solution is "consistent" with the core. Used internally by gapFillFastLP.
%
% Based on FASTCORE (Vlassis et al. 2014, PLoS Comput Biol 10:e1003424).
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
%     reaction carries flux in the L1-minimal solution. Empty if the LP is
%     infeasible (core reactions cannot all carry flux simultaneously).
%
% See also
% --------
% gapFillSwiftCore, gapFillFastLP

if nargin < 3 || isempty(epsilon)
    epsilon = 1e-4;
end

% ---- Resolve core reaction indices into original model ----
% coreOrigIdx: numeric indices into model.rxns
if iscell(coreRxns)
    coreOrigIdx = find(getIndexes(model, coreRxns, 'rxns', true));
elseif islogical(coreRxns)
    coreOrigIdx = find(coreRxns);
else
    coreOrigIdx = coreRxns(:);
end

% ---- Convert to irreversible form ----
% All fluxes become non-negative → L1 norm = sum of fluxes.
[irrevModel, ~, ~, irrev2rev] = convertToIrrev(model);
nIrrev = numel(irrevModel.rxns);
nMets  = numel(irrevModel.mets);

% Map core reactions to their irreversible equivalents.
% irrev2rev(k) gives the original-model index for irrev reaction k.
coreIdxIrrev = ismember(irrev2rev, coreOrigIdx);

% ---- Formulate LP ----
% Minimize sum of non-core fluxes (all fluxes >= 0 in irreversible form,
% so this is the L1 norm over non-core reactions).
%
%   min   sum_{j not in C} v_j
%   s.t.  S_irrev * v = 0
%         v_j >= epsilon   for j in C
%         lb <= v <= ub
%
lb = irrevModel.lb;
lb(coreIdxIrrev) = max(lb(coreIdxIrrev), epsilon);

c = ones(nIrrev, 1);
c(coreIdxIrrev) = 0;   % core reactions are not in the objective

prob.A      = irrevModel.S;
prob.a      = prob.A;
prob.b      = zeros(nMets, 1);
prob.csense = repmat('E', 1, nMets);
prob.c      = c;
prob.osense = 1;          % minimise
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
