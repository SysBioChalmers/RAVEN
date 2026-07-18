function [outModel, placement, addedTransports, exitFlag, report] = assignCompartments(model, GSS, reactionsToRelocate, varargin)
% assignCompartments  Assign reactions to compartments, certified by growth.
%
% Faithful MATLAB port of raven_toolbox.localization.assign_compartments. Places
% reactions into subcellular compartments from soft gene x compartment
% localization scores while keeping the result functional, WITHOUT ever putting
% a flux model inside the placement optimisation:
%
%   1. a flux-free placement master (maximise localization score,
%      mono-localisation): with no flux variable and no growth constraint, a
%      tolerance-rounded binary has nothing to leak into;
%   2. a structural confinement repair: reactions sharing a non-transportable
%      metabolite are co-located (or pinned to a fixed reaction's compartment),
%      and transportable pools a placement splits get star-topology transports;
%   3. materialised-FBA certification: the placement is confirmed functional by a
%      real solveLP on the model i_applyAssignment builds, over the primary
%      medium and every growth condition;
%   4. feedback on real failure: for a genuine placement gap, tighten the
%      placement and re-solve until certified or the round budget is hit.
%
% The model is NOT merged: each reaction keeps its current compartment, movable
% reactions are re-placed, and pinned reactions stay where they are (merge the
% model first, as a caller, to place onto a single-compartment draft).
%
% Parameters
% ----------
% model : struct
%     a RAVEN model with an objective (model.c) and grRules.
% GSS : struct
%     gene scoring structure (genes, compartments, scores) as from parseScores.
% reactionsToRelocate : cell
%     reaction ids to (re)place. Boundary, multi-compartment and objective
%     reactions are always pinned.
%
% Name-Value Arguments
% --------------------
% defaultCompartment : char
%     compartment transports route through (usually cytosol); must be in the
%     union of model and GSS compartments.
% multiCompartmentPenalty : double (default 0.5)
%     score cost per compartment a gene ends up in.
% minGrowth : double (default [] = 10%% of the unconstrained optimum)
%     required objective flux for certification.
% transportable : cell (default [] = all movable base metabolites)
%     metabolite NAMES that may receive transports.
% growthConditions : struct array (default [])
%     extra media to certify on, each with fields name, medium
%     (struct exchangeRxnId -> max uptake) and minGrowth.
% maxRounds : double (default 8)
%     budget on placement-tightening rounds.
% pruneTransports : logical (default true)
%     drop transports blocked in every certification medium.
% minimizeTransports : logical (default false)
%     after certifying, prune transports to those carrying flux in a
%     parsimonious-FBA solution, re-certifying (falls back if a needed
%     transport is pruned).
% biomassReaction : char (default [] = the largest objective-coefficient
%     reaction)
%     the reaction whose flux is the growth objective.
% baseMetabolite : char (default 'name')
%     the compartment-agnostic metabolite key: 'name' (metNames, the RAVEN
%     convention) or 'id' (model.mets).
% universal : struct (default [])
%     a template model; if a placement fails to certify because the primary
%     medium falls short, gaps are filled from it (via fillGaps).
% multiLocalize : logical (default false)
%     after certifying, add multi-compartment placements: for each placed
%     reaction, propose a second compartment its genes score highest for,
%     materialise the duplicate, and keep it only if a loopless FVA
%     (looplessFVA) shows it can carry flux in a biomass-supporting solution.
% multiLocalizeThreshold : double (default 0.7)
%     the gene score a second compartment needs for a duplicate to be proposed.
% multiLocalizeEps : double (default 1e-6)
%     the loopless flux a proposed duplicate must reach to be kept.
% transportCost : double (default 0.5)
%     accepted for signature compatibility but NOT used (placement is
%     flux-free and score-only).
% verbose : logical (default true)
%
% Returns
% -------
% outModel : struct
%     the compartmentalised model.
% placement : struct
%     .rxns and .compartment (cell of compartment ids per placed reaction).
% addedTransports : struct
%     .mets (base metabolite names) and .compartment for the added transports.
% exitFlag : double
%     1 = certified on every medium, -1 = could not place or did not certify.
% report : struct
%     .certified, .status, .growths (struct medium -> flux), .unplaced.
%
% See also
% --------
% predictLocalization, parseScores, gapFillMILP

p = parseRAVENargs(varargin, {'defaultCompartment',[]; 'transportCost',0.5; ...
    'multiCompartmentPenalty',0.5; 'minGrowth',[]; 'transportable',[]; ...
    'growthConditions',[]; 'maxRounds',8; 'pruneTransports',true; ...
    'minimizeTransports',false; 'biomassReaction',[]; 'baseMetabolite',[]; ...
    'universal',[]; 'multiLocalize',false; 'multiLocalizeThreshold',0.7; ...
    'multiLocalizeEps',1e-6; 'verbose',true});
defaultCompartment = char(p.defaultCompartment);
multiPen = p.multiCompartmentPenalty;
minGrowth = p.minGrowth;
growthConditions = p.growthConditions;
maxRounds = p.maxRounds;
pruneTransports = p.pruneTransports;
minimizeTransports = p.minimizeTransports;
universal = p.universal;
multiLocalize = p.multiLocalize;
mlThreshold = p.multiLocalizeThreshold;
mlEps = p.multiLocalizeEps;
verbose = p.verbose;

outModel = model;
placement = struct('rxns',{{}},'compartment',{{}});
addedTransports = struct('mets',{{}},'compartment',{{}});
exitFlag = -1;
report = struct('certified',false,'status','not_solved','growths',struct(),'unplaced',{{}});

if all(model.c == 0)
    error('RAVEN:badInput','model has no objective (set model.c).');
end
% compartments = union of model and score compartments
comps = unique([model.comps(:); GSS.compartments(:)], 'stable');
comps = sort(comps);
if ~ismember(defaultCompartment, comps)
    error('RAVEN:badInput','defaultCompartment ''%s'' not in the model/score compartments.', defaultCompartment);
end

% ---- scope: biomass, growth floor, movable/pinned, genes, transportable ----
sc = i_prepareScope(model, GSS, reactionsToRelocate, comps, minGrowth, p.transportable, p.biomassReaction, p.baseMetabolite);
minGrowth = sc.minGrowth;
report.unplaced = sc.unplaced;

if verbose
    fprintf('assignCompartments: %d movable, %d pinned reactions, %d genes, %d compartments.\n', ...
        numel(sc.movIdx), numel(sc.pinIdx), numel(sc.geneIdx), numel(comps));
end

% ---- place -> repair -> certify (-> tighten) loop ----
forced = containers.Map('KeyType','double','ValueType','double');  % movable local idx -> comp idx
groups = {};
gapPinned = containers.Map('KeyType','double','ValueType','logical');
seen = {};
best = struct('placeIdx',[],'trBase',[],'trComp',[],'growths',struct('primary',-1), ...
    'certified',false,'status','uncertified');

for round = 1:maxRounds
    sig = i_signature(forced, groups);
    if any(strcmp(seen, sig))
        break;  % configuration already tried: diagnosers are cycling
    end
    seen{end+1} = sig; %#ok<AGROW>

    placeIdx = i_placementMaster(model, sc, comps, multiPen, forced, groups, verbose);
    if isempty(placeIdx)
        report.status = 'infeasible';
        return;
    end

    % confinement fixpoint
    [newForced, newGroups, relaxed] = i_diagnoseConfinement(model, sc, comps, placeIdx, gapPinned);
    changed = false;
    fk = keys(newForced);
    for i = 1:numel(fk)
        k = fk{i};
        if (~isKey(forced,k) || forced(k) ~= newForced(k)) && ~isKey(gapPinned,k)
            forced(k) = newForced(k); changed = true;
        end
    end
    for gi = 1:numel(newGroups)
        if ~i_hasGroup(groups, newGroups{gi})
            groups{end+1} = newGroups{gi}; changed = true; %#ok<AGROW>
        end
    end
    if changed
        continue;
    end

    % transports + materialise + certify
    [trBase, trComp] = i_splitTransports(model, sc, comps, defaultCompartment, placeIdx, relaxed);
    outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trBase, trComp);
    [ok, growths] = i_certify(outModel, sc.biomassId, minGrowth, growthConditions);

    % feedback: gap-fill from the universal model if the primary medium fell short
    if ~ok && ~isempty(universal) && growths.primary < minGrowth - 1e-9
        gf = i_gapfill(outModel, universal, sc.biomassId, minGrowth);
        if ~isempty(gf.rxns)
            outModel = gf; [ok, growths] = i_certify(outModel, sc.biomassId, minGrowth, growthConditions);
        end
    end

    if ok
        if pruneTransports && ~isempty(trBase)
            keep = i_usableTransports(outModel, trComp, comps, sc.biomassId, growthConditions);
            trBase = trBase(keep); trComp = trComp(keep);
            outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trBase, trComp);
        end
        if minimizeTransports && ~isempty(trBase)
            [trBase, trComp, outModel] = i_minimizeTransports(model, sc, comps, ...
                defaultCompartment, placeIdx, trBase, trComp, minGrowth, growthConditions);
        end
        multiLoc = {};
        if multiLocalize
            [outModel, multiLoc] = i_enrichMultiloc(model, sc, comps, defaultCompartment, ...
                placeIdx, outModel, minGrowth, mlThreshold, mlEps, growthConditions);
        end
        exitFlag = 1; report.status = 'certified'; report.certified = true; report.growths = growths;
        report.multiLocalized = multiLoc;
        placement.rxns = model.rxns(sc.movIdx); placement.compartment = comps(placeIdx);
        addedTransports.mets = sc.baseNames(trBase); addedTransports.compartment = comps(trComp);
        return;
    end

    % keep best partial (largest primary growth)
    if growths.primary > best.growths.primary
        best = struct('placeIdx',placeIdx,'trBase',trBase,'trComp',trComp,'growths',growths, ...
            'certified',false,'status','uncertified');
    end

    % growth-gap feedback: pin the sole producer of a stranded biomass precursor
    fb = i_diagnoseGrowthGap(outModel, model, sc, comps, placeIdx);
    if ~isempty(fb) && (~isKey(forced,fb(1)) || forced(fb(1)) ~= fb(2))
        forced(fb(1)) = fb(2); gapPinned(fb(1)) = true;
        continue;
    end
    break;  % no tightening available
end

% honest uncertified result: return the best partial found
if ~isempty(best.placeIdx)
    outModel = i_applyAssignment(model, sc, comps, defaultCompartment, best.placeIdx, best.trBase, best.trComp);
    placement.rxns = model.rxns(sc.movIdx); placement.compartment = comps(best.placeIdx);
    addedTransports.mets = sc.baseNames(best.trBase); addedTransports.compartment = comps(best.trComp);
    report.growths = best.growths;
end
report.status = 'uncertified';
if verbose
    fprintf('assignCompartments: did not certify (growth %.4g < %.4g).\n', best.growths.primary, minGrowth);
end
end

% ============================================================ scope
function sc = i_prepareScope(model, GSS, relocate, comps, minGrowth, transportable, biomassReaction, baseMetabolite)
% biomass reaction: the largest-objective-coefficient reaction (or a named one)
if isempty(biomassReaction)
    [~, biomassIdx] = max(model.c);
else
    biomassIdx = find(strcmp(model.rxns, biomassReaction), 1);
    if isempty(biomassIdx)
        error('RAVEN:badInput','biomassReaction ''%s'' not in model.', biomassReaction);
    end
end
sc.biomassId = model.rxns{biomassIdx};

if isempty(minGrowth)
    sol = solveLP(model);
    if isempty(sol.f) || abs(sol.f) <= 0
        error('RAVEN:badInput','the draft model does not grow; pass minGrowth.');
    end
    minGrowth = 0.1 * abs(sol.f);
end
sc.minGrowth = minGrowth;

% base metabolite key (identifies the same species across compartments).
% Default 'name': RAVEN gives a metabolite a different id per compartment but a
% shared metName, so the name is the compartment-agnostic key (the equivalent
% of the Python default's compartment-suffix strip on cobra ids).
if isempty(baseMetabolite); baseMetabolite = 'name'; end
if strcmp(baseMetabolite,'name'); sc.base = model.metNames; else; sc.base = model.mets; end

% each reaction's single compartment ('' if boundary/multi-compartment)
sc.rxnComp = i_reactionCompartments(model);

toRelocate = ismember(model.rxns, relocate);
isBoundary = (sum(model.S ~= 0, 1)' == 1);
movable = toRelocate & ~isBoundary & ~cellfun(@isempty, sc.rxnComp);
movable(biomassIdx) = false;
% movable in sorted-reaction-id order (matching the Python master's variable
% order, so an equal-cost placement is broken the same way).
sc.movIdx = find(movable);
[~, ord] = sort(model.rxns(sc.movIdx));
sc.movIdx = sc.movIdx(ord);
sc.pinIdx = find(~movable);

% genes on movable reactions that are scored
[gInGSS, gssRow] = ismember(model.genes, GSS.genes);
geneOnMov = any(model.rxnGeneMat(sc.movIdx, :) ~= 0, 1)';
scope = gInGSS & geneOnMov;
sc.geneIdx = find(scope);
% genes in sorted-id order, so the placement MILP is built in a canonical
% order (identical to the Python master) and the solver sees the same problem.
[~, go] = sort(model.genes(sc.geneIdx));
sc.geneIdx = sc.geneIdx(go);
sc.score = zeros(numel(sc.geneIdx), numel(comps));
% align score columns (GSS.compartments) to comps
[~, colOf] = ismember(GSS.compartments, comps);
for gi = 1:numel(sc.geneIdx)
    row = GSS.scores(gssRow(sc.geneIdx(gi)), :);
    sc.score(gi, colOf(colOf>0)) = row(colOf>0);
end

% unplaced: movable reactions with genes but none scored
sc.unplaced = {};
for k = 1:numel(sc.movIdx)
    g = find(model.rxnGeneMat(sc.movIdx(k),:) ~= 0);
    if ~isempty(g) && ~any(ismember(g, sc.geneIdx))
        sc.unplaced{end+1,1} = model.rxns{sc.movIdx(k)}; %#ok<AGROW>
    end
end

% base metabolites: the same species across compartments shares a base key, so
% confinement and transports are keyed by base (not by per-compartment row),
% matching Python. baseId maps each metabolite row to its base index.
sc.metNames = model.metNames;
[sc.baseNames, baseRep, sc.baseId] = unique(sc.base, 'stable');
sc.nBase = numel(sc.baseNames);
sc.baseRep = baseRep;                                          % a met row per base
movBase = false(sc.nBase,1);
mrows = find(any(model.S(:, sc.movIdx) ~= 0, 2));
movBase(sc.baseId(mrows)) = true;                             % bases touched by movable rxns
if isnumeric(transportable) && isempty(transportable)
    sc.baseTransp = movBase;                                   % default: all movable bases
else
    sc.baseTransp = movBase & ismember(sc.baseNames, transportable);
end
end

% ============================================================ placement MILP
function placeIdx = i_placementMaster(model, sc, comps, multiPen, forced, groups, verbose)
% Flux-free score MILP. Variables x[mov,c] and y[gene,c], both binary.
% Objective: max sum score*y - multiPen*sum y (no reward on x).
movIdx = sc.movIdx; geneIdx = sc.geneIdx; score = sc.score;
nMov = numel(movIdx); nGene = numel(geneIdx); nC = numel(comps);
oX = 0;      nX = nMov*nC;
oY = oX+nX;  nY = nGene*nC;
nVar = oY+nY;
xCol = @(mi,ci) oX + (mi-1)*nC + ci;
yCol = @(gi,ci) oY + (gi-1)*nC + ci;

% Constraint rows are emitted in the exact order the Python master builds them,
% because the row order (like the column order) selects which of the degenerate
% co-optimal placements the solver returns -- a permutation of the rows changes
% the vertex. Order: per movable, its place row then its couple rows (in-scope
% genes on the reaction, in sorted-id order, x compartments); then per gene, its
% gene1 row then its has rows (x compartments); then the forced and colocation
% rows. Assembled as one triplet list with a running row index.
aR=[];aC=[];aV=[]; bVec=[]; cs=''; row=0;
for k=1:nMov
    row=row+1;                                            % place: sum_c x[k,c] = 1
    for ci=1:nC; aR(end+1)=row; aC(end+1)=xCol(k,ci); aV(end+1)=1; end %#ok<AGROW>
    bVec(end+1)=1; cs(end+1)='E'; %#ok<AGROW>
    gOn=find(model.rxnGeneMat(movIdx(k),geneIdx)~=0);     % in-scope genes on k, sorted-id order
    for gi=gOn
        for ci=1:nC                                      % couple: x[k,c] - y[gi,c] <= 0
            row=row+1;
            aR(end+1)=row; aC(end+1)=xCol(k,ci);  aV(end+1)=1;  %#ok<AGROW>
            aR(end+1)=row; aC(end+1)=yCol(gi,ci); aV(end+1)=-1; %#ok<AGROW>
            bVec(end+1)=0; cs(end+1)='L'; %#ok<AGROW>
        end
    end
end
for gi=1:nGene
    row=row+1;                                            % gene1: sum_c y[gi,c] >= 1
    for ci=1:nC; aR(end+1)=row; aC(end+1)=yCol(gi,ci); aV(end+1)=1; end %#ok<AGROW>
    bVec(end+1)=1; cs(end+1)='G'; %#ok<AGROW>
    rOfG=find(model.rxnGeneMat(movIdx,geneIdx(gi))~=0)';  % movable reactions on gi, movable order
    for ci=1:nC                                          % has: y[gi,c] - sum_r x[r,c] <= 0
        row=row+1; aR(end+1)=row; aC(end+1)=yCol(gi,ci); aV(end+1)=1; %#ok<AGROW>
        for k=rOfG; aR(end+1)=row; aC(end+1)=xCol(k,ci); aV(end+1)=-1; end %#ok<AGROW>
        bVec(end+1)=0; cs(end+1)='L'; %#ok<AGROW>
    end
end
fkeys=keys(forced);
for i=1:numel(fkeys)
    k=fkeys{i};                                          % force: x[k,forced(k)] = 1
    row=row+1; aR(end+1)=row; aC(end+1)=xCol(k,forced(k)); aV(end+1)=1; %#ok<AGROW>
    bVec(end+1)=1; cs(end+1)='E'; %#ok<AGROW>
end
for gi=1:numel(groups)
    mem=groups{gi};
    for j=1:numel(mem)-1
        a=mem(j); b=mem(j+1);
        for ci=1:nC                                      % colo: x[a,c] - x[b,c] = 0
            row=row+1;
            aR(end+1)=row; aC(end+1)=xCol(a,ci); aV(end+1)=1;  %#ok<AGROW>
            aR(end+1)=row; aC(end+1)=xCol(b,ci); aV(end+1)=-1; %#ok<AGROW>
            bVec(end+1)=0; cs(end+1)='E'; %#ok<AGROW>
        end
    end
end

c = zeros(nVar,1);
for gi=1:nGene
    for ci=1:nC; c(yCol(gi,ci)) = score(gi,ci) - multiPen; end
end

prob.A = sparse(aR,aC,aV,row,nVar);
prob.a = prob.A;
prob.b = bVec(:);
prob.csense = cs;
prob.c = -c; prob.osense = 1;
prob.lb = zeros(nVar,1); prob.ub = ones(nVar,1);
prob.vartype = repmat('B',1,nVar);

params.intTol = 1e-9; params.TimeLimit = 1000;
% deterministic, solver-independent solve: single thread and fixed seed remove
% thread-count/seed nondeterminism, MIPGap 0 forces the exact optimum rather
% than an early heuristic stop. With a canonically ordered problem this makes
% the placement reproducible and identical to the Python master.
params.Threads = 1; params.Seed = 0; params.MIPGap = 0;
sol = optimizeProb(prob, params, verbose);
if ~checkSolution(sol)
    placeIdx = []; return;
end
xval = sol.full(oX+(1:nX));
placeIdx = zeros(nMov,1);
for k=1:nMov
    [~,ci] = max(xval((k-1)*nC + (1:nC)));
    placeIdx(k) = ci;
end
end

% ============================================================ confinement
function [forced, groups, relaxed] = i_diagnoseConfinement(model, sc, comps, placeIdx, gapPinned)
% Find non-transportable metabolites a placement splits across compartments;
% force to a pinned compartment, co-locate all-movable sets, or relax
% (transport anyway) when pinned into >=2 comps or shared with a protected pin.
forced = containers.Map('KeyType','double','ValueType','double');
groups = {};
relaxed = containers.Map('KeyType','double','ValueType','logical');

placedComp = i_placedCompartments(model, sc, comps, placeIdx);   % per reaction, comp idx or 0
movLocal = containers.Map('KeyType','double','ValueType','double'); % rxn idx -> movable local idx
for k=1:numel(sc.movIdx); movLocal(sc.movIdx(k)) = k; end

% used{base} = compartments (and touching reactions) a non-transportable base
% appears in, keyed by base index.
usedComp = cell(sc.nBase,1); usedRxn = cell(sc.nBase,1);
allRxn = [sc.movIdx; sc.pinIdx];
for r = allRxn'
    comp = placedComp(r);
    if comp == 0; continue; end     % multi-compartment reaction: bridges pools
    for b = unique(sc.baseId(model.S(:,r)~=0))'
        if sc.baseTransp(b); continue; end
        usedComp{b}(end+1) = comp; usedRxn{b}(end+1) = r;
    end
end

for mm=1:sc.nBase
    cs = unique(usedComp{mm});
    if numel(cs) <= 1; continue; end
    touching = usedRxn{mm};
    isMov = arrayfun(@(r) isKey(movLocal,r), touching);
    movers = touching(isMov);
    % pinned compartments = comps where a pinned (non-movable) reaction touches mm
    pinnedComps = unique(arrayfun(@(i) placedComp(touching(i)), find(~isMov)));
    anyProtected = any(arrayfun(@(r) isKey(movLocal,r) && isKey(gapPinned,movLocal(r)), touching));
    if numel(pinnedComps) > 1 || anyProtected
        relaxed(mm) = true;
    elseif ~isempty(pinnedComps)
        for r = movers; forced(movLocal(r)) = pinnedComps(1); end
    elseif ~isempty(movers)
        g = sort(unique(arrayfun(@(r) movLocal(r), movers)));
        groups{end+1} = g(:)'; %#ok<AGROW>
    end
end
end

% ============================================================ split transports
function [trBase, trComp] = i_splitTransports(model, sc, comps, defaultCompartment, placeIdx, relaxed)
% Returns (baseIndex, compIndex) transports for transportable/relaxed bases a
% placement splits across compartments (routed through defaultCompartment).
defC = find(strcmp(comps, defaultCompartment), 1);
placedComp = i_placedCompartments(model, sc, comps, placeIdx);
allRxn = [sc.movIdx; sc.pinIdx];
used = cell(sc.nBase,1);
for r = allRxn'
    comp = placedComp(r);
    if comp == 0; continue; end
    for b = unique(sc.baseId(model.S(:,r)~=0))'
        if sc.baseTransp(b) || isKey(relaxed,b)
            used{b}(end+1) = comp;
        end
    end
end
trBase=[]; trComp=[];
for b=1:sc.nBase
    cs = unique(used{b});
    if numel(cs) <= 1; continue; end
    for ci = sort(cs)
        if ci ~= defC; trBase(end+1)=b; trComp(end+1)=ci; end %#ok<AGROW>
    end
end
trBase=trBase(:); trComp=trComp(:);
end

% ============================================================ certification
function [ok, growths] = i_certify(outModel, biomassId, minGrowth, growthConditions)
tol = 1e-9;
growths = struct();
growths.primary = i_growOn(outModel, biomassId, []);
ok = growths.primary >= minGrowth - tol;
if ~isempty(growthConditions) && isstruct(growthConditions)
    for i=1:numel(growthConditions)
        gc = growthConditions(i);
        g = i_growOn(outModel, biomassId, gc.medium);
        growths.(matlab.lang.makeValidName(gc.name)) = g;
        ok = ok && (g >= gc.minGrowth - tol);
    end
end
end

function g = i_growOn(model, biomassId, medium)
if ~isempty(medium); model = i_applyMedium(model, medium); end
bIdx = find(strcmp(model.rxns, biomassId), 1);
model.c = zeros(numel(model.rxns),1); model.c(bIdx) = 1;
sol = solveLP(model);
if isempty(sol.f); g = 0; else; g = abs(sol.f); end
end

function model = i_applyMedium(model, medium)
[~, exchIdx] = getExchangeRxns(model);
model.lb(exchIdx(model.lb(exchIdx) < 0)) = 0;
if ~isempty(medium) && isstruct(medium)
    f = fieldnames(medium);
    for i=1:numel(f)
        r = find(strcmp(model.rxns, f{i}), 1);
        if ~isempty(r); model.lb(r) = -abs(medium.(f{i})); end
    end
end
end

% ============================================================ transport pruning
function keep = i_usableTransports(outModel, trComp, comps, biomassId, growthConditions)
% Logical mask of transports that can carry flux in at least one certification
% medium (sound: an unusable reaction's removal cannot change any of those
% FBAs).
nTr = numel(trComp);
trId = cell(nTr,1);
for i=1:nTr
    trId{i} = ['tr_' num2str(i-1) '_' comps{trComp(i)}];
end
media = {[]};
if ~isempty(growthConditions) && isstruct(growthConditions)
    for i=1:numel(growthConditions); media{end+1} = growthConditions(i).medium; end %#ok<AGROW>
end
keep = false(nTr,1);
for mi=1:numel(media)
    m = outModel;
    if ~isempty(media{mi}); m = i_applyMedium(m, media{mi}); end
    bIdx = find(strcmp(m.rxns, biomassId), 1);
    m.c = zeros(numel(m.rxns),1); m.c(bIdx) = 1;
    idx = zeros(nTr,1);
    for i=1:nTr; j=find(strcmp(m.rxns,trId{i}),1); if ~isempty(j); idx(i)=j; end; end
    valid = idx>0;
    fl = false(nTr,1);
    fl(valid) = haveFlux(m, 'rxns', idx(valid));
    keep = keep | fl;
end
end

% ============================================================ universal gap-fill
function outModel = i_gapfill(outModel, universal, biomassId, minGrowth)
% Fill gaps from the universal model to restore the growth floor, using
% RAVEN's fillGaps. The added reactions are validated by the caller's
% certification FBA. On any failure the model is returned unchanged.
outModel.rxns = outModel.rxns;   % ensure struct is returned even on failure
try
    bIdx = find(strcmp(outModel.rxns, biomassId), 1);
    m = outModel; m.lb(bIdx) = minGrowth;   % require growth while filling
    [~, ~, addedRxns, newModel] = evalc('fillGaps(m, {universal}, ''useModelConstraints'', true, ''minGrowth'', minGrowth, ''verbose'', false)');
    if ~isempty(addedRxns)
        newModel.lb(strcmp(newModel.rxns, biomassId)) = outModel.lb(bIdx);  % restore biomass bound
        outModel = newModel;
    end
catch
    % gap-fill infeasible or unavailable: leave the model unchanged
end
end

% ============================================================ multi-localisation
function [outModel, multiLoc] = i_enrichMultiloc(model, sc, comps, defaultCompartment, placeIdx, outModel, minGrowth, threshold, eps, growthConditions)
% Add sound multi-compartment placements to a certified mono model: propose,
% per placed reaction, a second compartment its genes score >= threshold for;
% materialise every candidate as a duplicate; keep only those a loopless FVA
% shows can carry flux >= eps at the growth floor. Re-certifies; on failure or
% no survivor, returns the mono model.
multiLoc = {};
cand = [];  % [movableLocalIdx, compIdx]
for k=1:numel(sc.movIdx)
    genes = find(model.rxnGeneMat(sc.movIdx(k),:) ~= 0);
    gl = arrayfun(@(g) find(sc.geneIdx==g,1), genes, 'UniformOutput', false);
    gl = [gl{:}];
    if isempty(gl); continue; end
    for ci=1:numel(comps)
        if ci == placeIdx(k); continue; end
        if max(sc.score(gl, ci)) >= threshold
            cand = [cand; k, ci]; %#ok<AGROW>
        end
    end
end
if isempty(cand); return; end

% materialise all candidate duplicates onto the certified model
trial = outModel; dupId = cell(size(cand,1),1);
for i=1:size(cand,1)
    r = sc.movIdx(cand(i,1)); comp = comps{cand(i,2)};
    [trial, dupId{i}] = i_duplicateReaction(trial, model, r, comp);
end
trial = i_reconcileFields(trial);

% loopless FVA over the duplicates at the growth floor; keep flux-carriers
present = find(~cellfun(@isempty, dupId));
if isempty(present); return; end
ids = dupId(present);
[lo, hi] = looplessFVA(trial, ids, minGrowth);
keep = present(max(abs(lo), abs(hi)) >= eps);
if isempty(keep); return; end

% rebuild the model keeping only surviving duplicates, re-certify
enriched = outModel;
for i=keep'
    r = sc.movIdx(cand(i,1)); comp = comps{cand(i,2)};
    enriched = i_duplicateReaction(enriched, model, r, comp);
    multiLoc{end+1,1} = {model.rxns{r}, comp}; %#ok<AGROW>
end
enriched = i_reconcileFields(enriched);
if i_certify(enriched, sc.biomassId, minGrowth, growthConditions)
    outModel = enriched;
else
    multiLoc = {};   % safety net: duplicates only add capability, keep mono
end
end

function [outModel, dupId] = i_duplicateReaction(outModel, model, r, comp)
dupId = [model.rxns{r} '_' comp];
if any(strcmp(outModel.rxns, dupId)); dupId = ''; return; end
outModel.rxns{end+1,1} = dupId;
outModel.S(:, end+1) = 0;
for mm = find(model.S(:,r) ~= 0)'
    [outModel, newMet] = i_metInComp(outModel, model, mm, comp);
    outModel.S(newMet, end) = model.S(mm, r);
end
outModel.lb(end+1,1) = model.lb(r); outModel.ub(end+1,1) = model.ub(r);
if isfield(outModel,'rev'); outModel.rev(end+1,1) = model.rev(r); end
if isfield(outModel,'c'); outModel.c(end+1,1) = 0; end
if isfield(outModel,'rxnNames'); outModel.rxnNames{end+1,1} = dupId; end
if isfield(outModel,'grRules'); outModel.grRules{end+1,1} = model.grRules{r}; end
if isfield(outModel,'rxnGeneMat'); outModel.rxnGeneMat(end+1,:) = 0; end
end

% ============================================================ transport minimisation
function [trBase, trComp, outModel] = i_minimizeTransports(model, sc, comps, defaultCompartment, placeIdx, trBase, trComp, minGrowth, growthConditions)
% Prune the transports to those carrying flux in a parsimonious-FBA solution,
% then re-certify. Sound: a zero-flux reaction can be removed without changing
% the solution. Falls back to the input set if pFBA is infeasible or the
% pruned set regresses any medium.
outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trBase, trComp);
bIdx = find(strcmp(outModel.rxns, sc.biomassId), 1);
m = outModel; m.c = zeros(numel(m.rxns),1); m.c(bIdx) = 1;
sol = solveLP(m, 1);
if isempty(sol.x); return; end
carry = false(numel(trBase),1);
for i=1:numel(trBase)
    j = find(strcmp(outModel.rxns, ['tr_' num2str(i-1) '_' comps{trComp(i)}]), 1);
    if ~isempty(j) && abs(sol.x(j)) > 1e-7; carry(i) = true; end
end
newBase = trBase(carry); newComp = trComp(carry);
minModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, newBase, newComp);
ok = i_certify(minModel, sc.biomassId, minGrowth, growthConditions);
if ~ok; return; end                       % a needed transport was pruned: keep the safe set
trBase = newBase; trComp = newComp; outModel = minModel;
end

% ============================================================ growth-gap feedback
function fb = i_diagnoseGrowthGap(outModel, model, sc, comps, placeIdx)
% If a biomass precursor is blocked, pin its sole movable producer to biomass's
% compartment. Returns [movableLocalIdx, compIdx] or [].
fb = [];
bIdx = find(strcmp(outModel.rxns, sc.biomassId), 1);
if isempty(bIdx); return; end
bioComp = i_reactionCompartments(outModel); bioComp = bioComp{bIdx};
if isempty(bioComp); return; end
bcIdx = find(strcmp(comps, bioComp),1);
% precursors: base names biomass consumes
precRows = find(outModel.S(:,bIdx) < 0);
precNames = unique(outModel.metNames(precRows));
% blocked reactions in the materialised model
canFlux = haveFlux(outModel);
placedComp = i_placedCompartments(model, sc, comps, placeIdx);
for k=1:numel(sc.movIdx)
    r = sc.movIdx(k);
    ri = find(strcmp(outModel.rxns, model.rxns{r}),1);
    if isempty(ri) || canFlux(ri); continue; end
    touchNames = model.metNames(model.S(:,r)~=0);
    if any(ismember(touchNames, precNames)) && placedComp(r) ~= bcIdx
        fb = [k, bcIdx]; return;
    end
end
end

% ============================================================ materialisation
function outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trBase, trComp)
outModel = model;
for ci = 1:numel(comps)
    if ~ismember(comps{ci}, outModel.comps)
        outModel.comps{end+1,1} = comps{ci};
        if isfield(outModel,'compNames'); outModel.compNames{end+1,1} = comps{ci}; end
    end
end
for k = 1:numel(sc.movIdx)
    r = sc.movIdx(k); comp = comps{placeIdx(k)};
    mrows = find(model.S(:, r) ~= 0);
    for mm = mrows'
        [outModel, newMet] = i_metInComp(outModel, model, mm, comp);
        outModel.S(newMet, r) = model.S(mm, r);
        if ~isequal(newMet, mm); outModel.S(mm, r) = 0; end
    end
end
for k = 1:numel(trBase)
    outModel = i_addTransport(outModel, model, sc.baseRep(trBase(k)), comps{trComp(k)}, defaultCompartment, k-1);
end
outModel = i_reconcileFields(outModel);
end

function model = i_reconcileFields(model)
% Pad every registered rxn/met/gene/comp-indexed field to the current entity
% count, so the materialised model stays consistent for downstream functions
% (haveFlux, removeReactions, ...). Only the core fields are grown as
% reactions/metabolites are added; the rest are reconciled here.
reg = ravenModelFields();
counts = struct('rxn',numel(model.rxns),'met',numel(model.mets), ...
    'gene',numel(model.genes),'comp',numel(model.comps));
for i=1:numel(reg)
    f = reg(i).name;
    if ~isfield(model,f) || strcmp(f,'rxns') || strcmp(f,'mets') || ...
            strcmp(f,'genes') || strcmp(f,'comps'); continue; end
    n = counts.(reg(i).type);
    cur = size(model.(f),1);
    if cur >= n; continue; end
    def = reg(i).default;
    for j = cur+1:n
        if iscell(model.(f)); model.(f){j,1} = def; else; model.(f)(j,1) = def; end
    end
end
end

function [outModel, idx] = i_metInComp(outModel, model, srcMetRow, comp)
ci = find(strcmp(outModel.comps, comp), 1);
name = model.metNames{srcMetRow};
cand = find(strcmp(outModel.metNames, name) & outModel.metComps == ci, 1);
if ~isempty(cand); idx = cand; return; end
newId = [model.mets{srcMetRow} '__' comp];
outModel.mets{end+1,1} = newId;
outModel.metNames{end+1,1} = name;
outModel.metComps(end+1,1) = ci;
outModel.S(end+1, :) = 0;
if isfield(outModel,'b'); outModel.b(end+1,1) = 0; end
idx = numel(outModel.mets);
end

function outModel = i_addTransport(outModel, model, srcRow, comp, defaultCompartment, n)
[outModel, dRow] = i_metInComp(outModel, model, srcRow, defaultCompartment);
[outModel, cRow] = i_metInComp(outModel, model, srcRow, comp);
trId = ['tr_' num2str(n) '_' comp];
if any(strcmp(outModel.rxns, trId)); return; end
outModel.rxns{end+1,1} = trId;
outModel.S(:, end+1) = 0; outModel.S(dRow, end) = -1; outModel.S(cRow, end) = 1;
outModel.lb(end+1,1) = -1000; outModel.ub(end+1,1) = 1000;
if isfield(outModel,'rev'); outModel.rev(end+1,1) = 1; end
if isfield(outModel,'c'); outModel.c(end+1,1) = 0; end
if isfield(outModel,'rxnNames'); outModel.rxnNames{end+1,1} = trId; end
if isfield(outModel,'grRules'); outModel.grRules{end+1,1} = ''; end
if isfield(outModel,'rxnGeneMat'); outModel.rxnGeneMat(end+1,:) = 0; end
end

% ============================================================ helpers
function comp = i_reactionCompartments(model)
% Per reaction: the single compartment id of its metabolites, '' if boundary
% or multi-compartment.
comp = cell(numel(model.rxns),1);
for r=1:numel(model.rxns)
    cc = unique(model.metComps(model.S(:,r)~=0));
    if numel(cc)==1; comp{r} = model.comps{cc}; else; comp{r} = ''; end
end
end

function pc = i_placedCompartments(model, sc, comps, placeIdx)
% Per reaction index: compartment index. Movable -> placeIdx; pinned -> its own
% compartment; 0 if boundary/multi-compartment.
pc = zeros(numel(model.rxns),1);
for k=1:numel(sc.movIdx); pc(sc.movIdx(k)) = placeIdx(k); end
for r = sc.pinIdx'
    if ~isempty(sc.rxnComp{r})
        ci = find(strcmp(comps, sc.rxnComp{r}),1);
        if ~isempty(ci); pc(r) = ci; end
    end
end
end

function tf = i_hasGroup(groups, g)
tf = false; g = sort(g(:)');
for i=1:numel(groups)
    if isequal(sort(groups{i}(:)'), g); tf = true; return; end
end
end

function s = i_signature(forced, groups)
fk = sort(cell2mat(keys(forced)));
fs = '';
for i=1:numel(fk); fs = [fs sprintf('%d:%d,', fk(i), forced(fk(i)))]; end %#ok<AGROW>
gs = '';
for i=1:numel(groups); gs = [gs '|' num2str(sort(groups{i}(:)'))]; end %#ok<AGROW>
s = [fs '#' gs];
end
