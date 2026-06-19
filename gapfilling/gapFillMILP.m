function [addedRxns, reversedRxns, newModel, exitFlag] = gapFillMILP(model, universalModel, varargin)
% gapFillMILP  Objective-based gap-filling via a global growth-floor MILP.
%
% Finds the minimum-cost set of modifications to make the model's objective
% (biomass) >= minGrowth. Supports two repair mechanisms from Kumar et al.
% 2007 (Bioinformatics 23:1626-1635):
%   1. Directionality reversal of existing draft reactions (binary y_rev)
%   2. Addition of reactions from a universal database (binary y_db)
%
% The universal model should contain all candidate reactions, including
% exchange and transport reactions if those repair classes are desired.
%
% This is a SINGLE-level MILP (not bilevel). It is computationally more
% expensive than gapFillFastLP but supports the directionality-reversal
% repair mechanism and directly optimises for objective feasibility.
%
% Parameters
% ----------
% model : struct
%     draft RAVEN model to be gap-filled.
% universalModel : struct
%     universal reaction database (RAVEN model).
%
% Optional name-value parameters
% --------------------------------
% 'minGrowth' (default [])
%     minimum objective value to achieve. If empty, runs FBA on the merged
%     model and uses 10% of the maximum as the floor.
% 'weights' (default [1 2 3 4])
%     cost weights [w_rev w_db w_exch w_trans]. w_rev is the cost per
%     reversed draft reaction; w_db is the cost per added universal
%     reaction. w_exch and w_trans are reserved for future use when
%     exchange and transport reactions are tracked separately.
% 'bigM' (default 1000)
%     Big-M constant for coupling constraints. Must be larger than any
%     expected flux magnitude.
% 'params' (default [])
%     solver parameter struct passed to optimizeProb. Fields such as
%     params.TimeLimit and params.intTol can be used to tune the MILP.
% 'verbose' (default true)
%     print progress messages.
%
% Returns
% -------
% addedRxns : cell
%     reaction IDs from universalModel that are selected by the MILP.
% reversedRxns : cell
%     draft reaction IDs whose directionality is reversed by the MILP.
% newModel : struct
%     model with addedRxns incorporated and reversedRxns bounds updated
%     (lb set to -bigM for reversed reactions).
% exitFlag : double
%     1 = optimal, -1 = infeasible, -2 = timeout or other failure.
%
% See also
% --------
% gapFillFastLP, gapFillTopological, fillGaps

% ---- Parse optional arguments ----
p = inputParser();
addParameter(p, 'minGrowth', [],  @(x) isempty(x) || isnumeric(x));
addParameter(p, 'weights',   [1 2 3 4], @isnumeric);
addParameter(p, 'bigM',      1000,      @isnumeric);
addParameter(p, 'params',    [],        @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose',   true,      @islogical);
parse(p, varargin{:});
minGrowth = p.Results.minGrowth;
weights   = p.Results.weights;
M         = p.Results.bigM;
params    = p.Results.params;
verbose   = p.Results.verbose;

w_rev = weights(1);
w_db  = weights(2);

% ---- Initialise outputs ----
addedRxns    = {};
reversedRxns = {};
newModel     = model;
exitFlag     = -1;

% Check that the model has an objective
if all(model.c == 0)
    error('gapFillMILP: model has no objective function. Set model.c.');
end

% ---- Merge draft + universal database ----
mergedModel = mergeModels({model; universalModel});
nMergedMets = numel(mergedModel.mets);
nMergedRxns = numel(mergedModel.rxns);

% Identify draft reactions in the merged model
isDraft = ismember(mergedModel.rxns, model.rxns);
isUniv  = ~isDraft;

% Universal reactions that can carry flux (non-trivially blocked)
canCarry = (mergedModel.ub > 0 | mergedModel.lb < 0);
isUniv   = isUniv & canCarry;
univIdx  = find(isUniv);     % positions in merged model
nUniv    = numel(univIdx);

% Reversal candidates: draft reactions that are currently irreversible
% (lb >= 0). Reversing their directionality means allowing negative flux.
revCandMask = isDraft & mergedModel.lb >= 0;
revCandIdx  = find(revCandMask);   % positions in merged model
nRev        = numel(revCandIdx);

if verbose
    fprintf('gapFillMILP: merged model has %d mets, %d rxns (%d draft, %d universal).\n', ...
        nMergedMets, nMergedRxns, sum(isDraft), nUniv);
    fprintf('gapFillMILP: %d reversal candidates, %d database candidates.\n', nRev, nUniv);
end

% ---- Determine minGrowth ----
if isempty(minGrowth)
    % Run FBA on merged model (unconstrained) to get max growth.
    probFBA.A      = mergedModel.S;
    probFBA.a      = probFBA.A;
    probFBA.b      = zeros(nMergedMets, 1);
    probFBA.csense = repmat('E', 1, nMergedMets);
    probFBA.c      = -mergedModel.c;   % negate: osense=1 minimises
    probFBA.osense = 1;
    probFBA.lb     = mergedModel.lb;
    probFBA.ub     = mergedModel.ub;
    probFBA.vartype = repmat('C', 1, nMergedRxns);
    solFBA = optimizeProb(probFBA, params, false);
    if ~checkSolution(solFBA)
        if verbose
            fprintf('gapFillMILP: FBA on merged model is infeasible — cannot gap-fill.\n');
        end
        return;
    end
    maxGrowth = mergedModel.c' * solFBA.full(1:nMergedRxns);
    if maxGrowth <= 0
        if verbose
            fprintf('gapFillMILP: merged model objective is non-positive (%.4g). Cannot gap-fill.\n', ...
                maxGrowth);
        end
        return;
    end
    minGrowth = 0.1 * maxGrowth;
    if verbose
        fprintf('gapFillMILP: setting minGrowth = %.4g (10%% of max %.4g).\n', minGrowth, maxGrowth);
    end
end

% ---- Build MILP ----
% Variable ordering:
%   cols 1 : nMergedRxns          → v (all fluxes, continuous)
%   cols nMergedRxns+1 : +nRev    → y_rev (binary, directionality reversal)
%   cols nMergedRxns+nRev+1 : end → y_db  (binary, universal reactions)
nVar = nMergedRxns + nRev + nUniv;

% --- Variable bounds ---
lb = mergedModel.lb;
ub = mergedModel.ub;

% Relax reversal candidates to allow negative flux (Big-M on lower bound)
lb(revCandIdx) = -M;

% Universal reactions: allow them to be forced to 0 (extend lb toward 0)
lb(univIdx) = min(mergedModel.lb(univIdx), 0);

% Binary variable bounds
lb = [lb; zeros(nRev + nUniv, 1)];
ub = [ub; ones(nRev + nUniv, 1)];

% --- Objective ---
c = zeros(nVar, 1);
c(nMergedRxns + 1 : nMergedRxns + nRev)  = w_rev;  % y_rev cost
c(nMergedRxns + nRev + 1 : end)           = w_db;   % y_db cost

% --- Constraint blocks ---
% (1) Steady state: S * v = 0   [nMergedMets rows, 'E']
A_SS = [mergedModel.S, sparse(nMergedMets, nRev + nUniv)];

% (2) Growth floor: c' * v >= minGrowth   [1 row, 'G']
A_growth = [sparse(mergedModel.c'), sparse(1, nRev + nUniv)];

% (3) Universal upper coupling: v_k - ub_k * y_db_k <= 0   [nUniv rows, 'L']
%     (forces v_k <= 0 when y_db_k = 0)
yDbCols = nMergedRxns + nRev + (1:nUniv);
A_uniUp = sparse(...
    [1:nUniv, 1:nUniv], ...
    [univIdx(:)', yDbCols], ...
    [ones(1, nUniv), -mergedModel.ub(univIdx)'], ...
    nUniv, nVar);

% (4) Universal lower coupling: v_k - lb_k * y_db_k >= 0   [nUniv rows, 'G']
%     (forces v_k >= 0 when y_db_k = 0; together with (3): v_k = 0 when y_db_k = 0)
A_uniLo = sparse(...
    [1:nUniv, 1:nUniv], ...
    [univIdx(:)', yDbCols], ...
    [ones(1, nUniv), -mergedModel.lb(univIdx)'], ...
    nUniv, nVar);

% (5) Reversal coupling: v_j + M * y_rev_j >= 0   [nRev rows, 'G']
%     (forces v_j >= 0 when y_rev_j = 0, i.e., preserves original irreversibility)
yRevCols = nMergedRxns + (1:nRev);
A_rev = sparse(...
    [1:nRev, 1:nRev], ...
    [revCandIdx(:)', yRevCols], ...
    [ones(1, nRev), M * ones(1, nRev)], ...
    nRev, nVar);

% --- Assemble ---
prob.A = [A_SS; A_growth; A_uniUp; A_uniLo; A_rev];
prob.a = prob.A;
prob.b = [zeros(nMergedMets, 1); minGrowth; ...
          zeros(nUniv, 1); zeros(nUniv, 1); zeros(nRev, 1)];
prob.csense = [repmat('E', 1, nMergedMets), 'G', ...
               repmat('L', 1, nUniv), ...
               repmat('G', 1, nUniv + nRev)];
prob.c      = c;
prob.osense = 1;   % minimise weighted sum of binary variables
prob.lb     = lb;
prob.ub     = ub;
prob.vartype = repmat('C', 1, nVar);
prob.vartype(nMergedRxns + 1 : end) = 'B';   % y_rev and y_db are binary

if isempty(params)
    params.intTol   = 1e-9;
    params.TimeLimit = 600;
end

% ---- Solve ----
if verbose
    fprintf('gapFillMILP: solving MILP (%d vars, %d binary, %d constraints)...\n', ...
        nVar, nRev + nUniv, size(prob.A, 1));
end
sol = optimizeProb(prob, params, verbose);

if ~checkSolution(sol)
    if verbose
        fprintf('gapFillMILP: MILP infeasible or timed out.\n');
    end
    exitFlag = -1;
    return;
end

% ---- Extract solution ----
binThreshold = 1e-3;
yRev = sol.full(nMergedRxns + 1 : nMergedRxns + nRev);
yDb  = sol.full(nMergedRxns + nRev + 1 : end);

revSelected = revCandIdx(yRev > binThreshold);
dbSelected  = univIdx(yDb > binThreshold);

% Map back to original IDs
reversedRxns = mergedModel.rxns(revSelected);
addedRxns    = mergedModel.rxns(dbSelected);
% Keep only universal reaction IDs (filtered by name)
addedRxns = addedRxns(ismember(addedRxns, universalModel.rxns));

if verbose
    fprintf('gapFillMILP: reversed %d draft reaction(s), added %d universal reaction(s).\n', ...
        numel(reversedRxns), numel(addedRxns));
end

% ---- Build new model ----
newModel = model;

% Apply directionality reversals: allow negative flux
for k = 1:numel(reversedRxns)
    rxnPos = find(getIndexes(newModel, reversedRxns(k), 'rxns', true));
    if ~isempty(rxnPos)
        newModel.lb(rxnPos) = -M;
        newModel.rev(rxnPos) = 1;
    end
end

% Add universal reactions
if ~isempty(addedRxns)
    trimmedUniv = removeReactions(universalModel, ...
        setdiff(universalModel.rxns, addedRxns));
    newModel = mergeModels({newModel; trimmedUniv});
end

exitFlag = 1;
end
