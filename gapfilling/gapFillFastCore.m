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
[irrevModel, matchRev, ~, irrev2rev] = convertToIrrev(model);
nIrrev = numel(irrevModel.rxns);
nMets  = numel(irrevModel.mets);

% Map core reactions to their irreversible equivalents. convertToIrrev
% keeps the forward copy at the reaction's original index, so coreOrigIdx
% already are that copy's indices.
coreIdxIrrev = false(nIrrev, 1);
coreIdxIrrev(coreOrigIdx) = true;

% A reversible core reaction also has a reverse copy elsewhere in
% irrevModel, with identical but oppositely-signed stoichiometry. Forcing
% both copies to >= epsilon would let them cancel out in S_irrev*v=0
% (v_fwd=v_bwd=epsilon nets to zero for every metabolite the reaction
% touches), satisfying the "core must carry flux" requirement with a
% self-contained loop that needs no support from the rest of the network,
% regardless of whether the reaction has any real connectivity. Instead,
% require flux on one copy at a time: try forward first, and if that is
% infeasible, retry with reverse forced instead, so a reaction whose only
% genuine support is in the reverse direction is still recognised.
isRevCore   = matchRev(coreOrigIdx) > 0;
fwdCoreCopy = coreOrigIdx(isRevCore);          % forward copy of reversible core rxns
revCoreCopy = matchRev(coreOrigIdx(isRevCore)); % their reverse copy
irrCoreCopy = coreOrigIdx(~isRevCore);          % irreversible core rxns: single copy

c = ones(nIrrev, 1);
c(coreIdxIrrev) = 0;   % core reactions are not in the objective

prob.A      = irrevModel.S;
prob.a      = prob.A;
prob.b      = zeros(nMets, 1);
prob.csense = repmat('E', 1, nMets);
prob.c      = c;
prob.osense = 1;          % minimise
prob.vartype = repmat('C', 1, nIrrev);

% ---- Formulate LP ----
% Minimize sum of non-core fluxes (all fluxes >= 0 in irreversible form,
% so this is the L1 norm over non-core reactions).
%
%   min   sum_{j not in C} v_j
%   s.t.  S_irrev * v = 0
%         v_j >= epsilon   for j in C
%         lb <= v <= ub
%
forceCopy   = {fwdCoreCopy, revCoreCopy};
disableCopy = {revCoreCopy, fwdCoreCopy};
nPasses     = 1 + ~isempty(fwdCoreCopy);   % a second pass only differs when a reversible core rxn exists
for pass = 1:nPasses
    lb = irrevModel.lb;
    lb(irrCoreCopy)         = max(lb(irrCoreCopy), epsilon);
    lb(forceCopy{pass})     = max(lb(forceCopy{pass}), epsilon);
    ub = irrevModel.ub;
    ub(disableCopy{pass})   = 0;
    prob.lb = lb;
    prob.ub = ub;

    sol = optimizeProb(prob, [], false);
    if checkSolution(sol)
        break;
    end
end

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
