function [addedRxns, newModel, cannotConnect] = gapFillFastLP(model, universalModel, varargin)
% gapFillFastLP  LP-based gap-filling (fastGapFill / swiftGapFill).
%
% Merges the draft model with a universal reaction database, then for each
% blocked draft reaction calls an LP subroutine (FASTCORE or SWIFTCORE) to
% find the minimum set of database reactions that makes the blocked reaction
% carry flux. The union of all activated database reactions across all blocked
% reactions is returned as addedRxns.
%
% Based on fastGapFill (Thiele et al. 2014, PLoS Comput Biol 10:e1003515)
% and swiftGapFill (Tefagh & Boyd 2020, BMC Bioinformatics 21:23).
%
% Parameters
% ----------
% model : struct
%     draft RAVEN model to be gap-filled.
% universalModel : struct
%     universal reaction database (RAVEN model). The caller is responsible
%     for building this, e.g. via getKEGGModelForOrganism with broad scope.
%
% Optional name-value parameters
% --------------------------------
% 'epsilon' (default 1e-4)
%     minimum flux threshold for a reaction to be considered active.
% 'variant' (default 'fast')
%     'fast'  — use FASTCORE L1-norm LP (gapFillFastCore)
%     'swift' — use SWIFTCORE single-LP (gapFillSwiftCore); faster but
%               stochastic (results may vary between runs/solvers).
% 'verbose' (default true)
%     print progress messages.
%
% Returns
% -------
% addedRxns : cell
%     reaction IDs from universalModel that should be added to rescue
%     blocked reactions.
% newModel : struct
%     copy of model with addedRxns incorporated.
% cannotConnect : cell
%     blocked draft reaction IDs that remain blocked even after augmenting
%     with the full universalModel (cannot be rescued by the database).
%
% See also
% --------
% gapFillFastCore, gapFillSwiftCore, gapFillMILP, fillGaps

% ---- Parse optional arguments ----
p = inputParser();
addParameter(p, 'epsilon', 1e-4,   @isnumeric);
addParameter(p, 'variant', 'fast', @ischar);
addParameter(p, 'verbose', true,   @islogical);
parse(p, varargin{:});
epsilon = p.Results.epsilon;
variant = lower(p.Results.variant);
verbose = p.Results.verbose;

% ---- Initialise outputs ----
addedRxns    = {};
cannotConnect = {};
newModel     = model;

% ---- Identify blocked reactions in the draft model ----
blocked = find(~haveFlux(model, epsilon));
if isempty(blocked)
    if verbose
        fprintf('gapFillFastLP: all reactions carry flux; nothing to repair.\n');
    end
    return;
end
nBlocked = numel(blocked);
if verbose
    fprintf('gapFillFastLP: %d blocked reaction(s) found in draft model.\n', nBlocked);
end

% ---- Merge draft + universal database ----
mergedModel = mergeModels({model; universalModel});

% For each draft reaction, find its position in the merged model.
[inMerged, draftPosInMerged] = ismember(model.rxns, mergedModel.rxns);
if any(~inMerged)
    warning('gapFillFastLP: %d draft reaction(s) not found in merged model — skipped.', sum(~inMerged));
end

% ---- Determine which blocked reactions can carry flux in the merged model ----
mergedHasFlux = haveFlux(mergedModel, epsilon);
draftCanCarry = mergedHasFlux(draftPosInMerged);   % indexed by draft position

% Blocked reactions that still cannot carry flux even with the full DB.
blockedCanCarry = draftCanCarry(blocked);
cannotConnectIdx = blocked(~blockedCanCarry);
cannotConnect    = model.rxns(cannotConnectIdx);
rescuable        = blocked(blockedCanCarry);

if verbose && ~isempty(cannotConnect)
    fprintf('gapFillFastLP: %d reaction(s) cannot be rescued by the universal database.\n', ...
        numel(cannotConnect));
end

if isempty(rescuable)
    if verbose
        fprintf('gapFillFastLP: no rescuable blocked reactions; nothing to add.\n');
    end
    return;
end
if verbose
    fprintf('gapFillFastLP: %d/%d blocked reaction(s) are rescuable; running %s...\n', ...
        numel(rescuable), nBlocked, variant);
end

% ---- For each rescuable blocked reaction, run the LP subroutine ----
% Accumulate which merged-model reactions are activated.
isUniversal   = ~ismember(mergedModel.rxns, model.rxns);
candidateSet  = false(numel(mergedModel.rxns), 1);

nSkipped = 0;
for i = 1:numel(rescuable)
    coreIdxInMerged = draftPosInMerged(rescuable(i));
    if coreIdxInMerged == 0
        nSkipped = nSkipped + 1;
        continue;
    end

    if strcmp(variant, 'swift')
        active = gapFillSwiftCore(mergedModel, coreIdxInMerged, epsilon);
    else
        active = gapFillFastCore(mergedModel, coreIdxInMerged, epsilon);
    end

    if any(active)
        candidateSet = candidateSet | active;
    else
        nSkipped = nSkipped + 1;
        if verbose
            fprintf('gapFillFastLP: LP infeasible for %s; skipping.\n', model.rxns{rescuable(i)});
        end
    end
end

if verbose && nSkipped > 0
    fprintf('gapFillFastLP: %d reaction(s) skipped (infeasible LP).\n', nSkipped);
end

% ---- Extract universal reactions that appear in any solution ----
addedMask = candidateSet & isUniversal;
addedRxns = mergedModel.rxns(addedMask);

if isempty(addedRxns)
    if verbose
        fprintf('gapFillFastLP: no universal reactions needed (all gaps already resolved).\n');
    end
    return;
end

% ---- Build new model: draft + selected universal reactions ----
toRemove = find(isUniversal & ~addedMask);
if ~isempty(toRemove)
    newModel = removeReactions(mergedModel, mergedModel.rxns(toRemove));
else
    newModel = mergedModel;
end

if verbose
    fprintf('gapFillFastLP: added %d reaction(s) from universal database.\n', numel(addedRxns));
end
end
