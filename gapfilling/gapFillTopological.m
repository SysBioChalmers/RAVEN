function result = gapFillTopological(model, universalModel, varargin)
% gapFillTopological  BFS metabolite-producibility pre-screening.
%
% Identifies which metabolites in the draft model cannot be produced from a
% given set of seed metabolites (medium), using a graph-based fixed-point
% "scope" computation (no solver required). For each unreachable metabolite,
% finds candidate reactions in the universal database that produce it
% directly.
%
% This is a computationally cheap pre-screening step that prunes the
% universal database before a more expensive gap-filling solve. It is
% inspired by the topological analysis used in Meneco (Prigent et al. 2017,
% PLoS Comput Biol).
%
% Algorithm
% ---------
% Starting from seed metabolites (available in the medium), iteratively
% mark reactions as "fireable" when all their substrates are reachable,
% then mark their products as newly reachable. Continue until no new
% metabolite is reached (fixed point). Reversible reactions (lb < 0) are
% also considered in their reverse direction.
%
% Parameters
% ----------
% model : struct
%     draft RAVEN model.
% universalModel : struct
%     universal reaction database (RAVEN model).
%
% Optional name-value parameters
% --------------------------------
% 'seeds' (default [])
%     Metabolite IDs (cell array) available from the medium. Default:
%     metabolites involved in uptake exchange reactions (lb < 0).
% 'targets' (default [])
%     Metabolite IDs (cell array) that should be produced. Default:
%     substrates of the objective (biomass) reaction.
% 'verbose' (default true)
%     print a summary of the results.
%
% Returns
% -------
% result : struct with fields:
%   .reachableMets   — logical vector (length nMets) marking producible mets
%   .blockedMets     — cell array of metabolite IDs that cannot be reached
%   .candidateRxns   — n-by-1 cell array of universal rxn ID lists, one
%                      cell per blocked metabolite (same order as blockedMets)
%   .pruningFraction — fraction of universal reactions pruned (not needed
%                      for any blocked metabolite)
%
% See also
% --------
% gapFillFastLP, gapFillMILP, fillGaps

% ---- Parse optional arguments ----
p = inputParser();
addParameter(p, 'seeds',   [], @(x) isempty(x) || iscell(x));
addParameter(p, 'targets', [], @(x) isempty(x) || iscell(x));
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});
seeds   = p.Results.seeds;
targets = p.Results.targets;
verbose = p.Results.verbose;

nMets = numel(model.mets);
nRxns = numel(model.rxns);

% ---- Identify seed metabolites ----
if isempty(seeds)
    % Default: metabolites in exchange reactions that allow uptake (lb < 0)
    [~, exchIdx] = getExchangeRxns(model);
    uptakeIdx = exchIdx(model.lb(exchIdx) < 0);
    if isempty(uptakeIdx)
        warning(['gapFillTopological: no uptake exchange reactions found. ' ...
            'Provide seeds manually via the ''seeds'' option.']);
        seedMetIdx = [];
    else
        [seedRows, ~] = find(model.S(:, uptakeIdx) ~= 0);
        seedMetIdx = unique(seedRows);
    end
else
    [~, seedMetIdx] = ismember(seeds, model.mets);
    seedMetIdx = seedMetIdx(seedMetIdx > 0);
end

% ---- Identify target metabolites ----
if isempty(targets)
    % Default: substrates of the objective (biomass) reaction
    objRxns = find(model.c ~= 0);
    if isempty(objRxns)
        error(['gapFillTopological: no objective function defined. ' ...
            'Provide targets manually via the ''targets'' option.']);
    end
    biomassIdx = objRxns(1);
    targetMetIdx = find(model.S(:, biomassIdx) ~= 0);
else
    [~, targetMetIdx] = ismember(targets, model.mets);
    targetMetIdx = targetMetIdx(targetMetIdx > 0);
end

% ---- Build adjacency lists ----
% subOf{m}  = reaction indices where metabolite m is a substrate (S < 0)
% prodOf{m} = reaction indices where metabolite m is a product  (S > 0)
% rxnSubs{j} = substrate metabolite indices of reaction j
% rxnProds{j} = product metabolite indices of reaction j
%
% Use find() on the sparse S directly to avoid materialising a dense logical
% matrix: extract all nonzeros then partition by sign.
[allRows, allCols, allVals] = find(model.S);
subMask    = allVals < 0;
proMask    = allVals > 0;
subMetRows = allRows(subMask);
subRxnCols = allCols(subMask);
proMetRows = allRows(proMask);
proRxnCols = allCols(proMask);

subOf  = cell(nMets, 1);
prodOf = cell(nMets, 1);
for k = 1:numel(subMetRows)
    subOf{subMetRows(k)}(end+1) = subRxnCols(k);
end
for k = 1:numel(proMetRows)
    prodOf{proMetRows(k)}(end+1) = proRxnCols(k);
end

rxnSubs  = cell(nRxns, 1);
rxnProds = cell(nRxns, 1);
for k = 1:numel(subMetRows)
    rxnSubs{subRxnCols(k)}(end+1) = subMetRows(k);
end
for k = 1:numel(proMetRows)
    rxnProds{proRxnCols(k)}(end+1) = proMetRows(k);
end

isRev = model.lb < 0;

% ---- BFS scope computation ----
% Uses a countdown approach: a reaction fires when all its substrates are
% reachable. Each time a metabolite becomes reachable, decrement the
% substrate count of all reactions it participates in.
subCountFwd = cellfun(@numel, rxnSubs);   % remaining unreachable substrates (fwd)
subCountRev = cellfun(@numel, rxnProds);  % remaining unreachable products  (rev)

reachable = false(nMets, 1);
firedFwd  = false(nRxns, 1);
firedRev  = false(nRxns, 1);

% Initialize queue with seed metabolites
queue  = seedMetIdx(:)';
qHead  = 1;
reachable(seedMetIdx) = true;

while qHead <= numel(queue)
    m = queue(qHead);
    qHead = qHead + 1;

    % --- Forward direction ---
    % m is a substrate of rxn j (S(m,j) < 0)
    for j = subOf{m}
        subCountFwd(j) = subCountFwd(j) - 1;
        if subCountFwd(j) == 0 && ~firedFwd(j)
            firedFwd(j) = true;
            for p = rxnProds{j}
                if ~reachable(p)
                    reachable(p) = true;
                    queue(end+1) = p;
                end
            end
        end
    end

    % --- Reverse direction (reversible reactions only) ---
    % m is a product of rxn j (S(m,j) > 0), which in reverse means m is a
    % substrate; if all products are reachable, the reaction can fire in
    % reverse and produce the substrates.
    for j = prodOf{m}
        if isRev(j)
            subCountRev(j) = subCountRev(j) - 1;
            if subCountRev(j) == 0 && ~firedRev(j)
                firedRev(j) = true;
                for s = rxnSubs{j}
                    if ~reachable(s)
                        reachable(s) = true;
                        queue(end+1) = s;
                    end
                end
            end
        end
    end
end

% ---- Report blocked target metabolites ----
targetReachable = reachable(targetMetIdx);
blockedMetIdx   = targetMetIdx(~targetReachable);
blockedMets     = model.mets(blockedMetIdx);

if verbose
    fprintf('gapFillTopological: %d/%d target metabolites are reachable.\n', ...
        sum(targetReachable), numel(targetMetIdx));
    if ~isempty(blockedMets)
        fprintf('gapFillTopological: %d blocked target metabolite(s).\n', numel(blockedMets));
    end
end

% ---- Find candidate universal reactions for each blocked metabolite ----
candidateRxns  = cell(numel(blockedMets), 1);
allCandidates  = {};

for k = 1:numel(blockedMetIdx)
    m     = blockedMetIdx(k);
    metId = model.mets{m};

    % Find this metabolite in the universal model
    univMetIdx = find(strcmp(universalModel.mets, metId));
    if isempty(univMetIdx)
        candidateRxns{k} = {};
        continue;
    end
    univMetIdx = univMetIdx(1);   % take first match; IDs should be unique

    univSRow = full(universalModel.S(univMetIdx, :));

    % Forward production: rxns with S(m,j) > 0
    fwdProdMask = univSRow > 0;
    % Reverse production: rxns with S(m,j) < 0 that can carry negative flux
    revProdMask = univSRow < 0 & (universalModel.lb' < 0);

    candidateMask         = fwdProdMask | revProdMask;
    candidateRxns{k}      = universalModel.rxns(candidateMask);
    allCandidates         = [allCandidates; candidateRxns{k}];
end

% ---- Pruning fraction ----
uniqueCandidates   = unique(allCandidates);
pruningFraction    = 1 - numel(uniqueCandidates) / numel(universalModel.rxns);

if verbose
    fprintf('gapFillTopological: %d/%d universal reactions pruned (%.0f%%).\n', ...
        numel(universalModel.rxns) - numel(uniqueCandidates), ...
        numel(universalModel.rxns), pruningFraction * 100);
end

% ---- Pack output ----
result.reachableMets     = reachable;
result.blockedMets       = blockedMets;
result.candidateRxns     = candidateRxns;
result.pruningFraction   = pruningFraction;
end
