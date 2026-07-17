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
    'growthConditions',[]; 'maxRounds',8; 'pruneTransports',true; 'verbose',true});
defaultCompartment = char(p.defaultCompartment);
multiPen = p.multiCompartmentPenalty;
minGrowth = p.minGrowth;
growthConditions = p.growthConditions;
maxRounds = p.maxRounds;
pruneTransports = p.pruneTransports;
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
sc = i_prepareScope(model, GSS, reactionsToRelocate, comps, defaultCompartment, minGrowth, p.transportable, verbose);
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
best = struct('placeIdx',[],'trMet',[],'trComp',[],'growths',struct('primary',-1), ...
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
    [trMet, trComp] = i_splitTransports(model, sc, comps, defaultCompartment, placeIdx, relaxed);
    outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trMet, trComp);
    [ok, growths] = i_certify(outModel, sc.biomassId, minGrowth, growthConditions);

    if ok
        if pruneTransports && ~isempty(trMet)
            [trMet, trComp] = i_usableTransports(outModel, trMet, trComp, comps, sc.biomassId, growthConditions);
            outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trMet, trComp);
        end
        exitFlag = 1; report.status = 'certified'; report.certified = true; report.growths = growths;
        placement.rxns = model.rxns(sc.movIdx); placement.compartment = comps(placeIdx);
        addedTransports.mets = sc.metNames(trMet); addedTransports.compartment = comps(trComp);
        return;
    end

    % keep best partial (largest primary growth)
    if growths.primary > best.growths.primary
        best = struct('placeIdx',placeIdx,'trMet',trMet,'trComp',trComp,'growths',growths, ...
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
    outModel = i_applyAssignment(model, sc, comps, defaultCompartment, best.placeIdx, best.trMet, best.trComp);
    placement.rxns = model.rxns(sc.movIdx); placement.compartment = comps(best.placeIdx);
    addedTransports.mets = sc.metNames(best.trMet); addedTransports.compartment = comps(best.trComp);
    report.growths = best.growths;
end
report.status = 'uncertified';
if verbose
    fprintf('assignCompartments: did not certify (growth %.4g < %.4g).\n', best.growths.primary, minGrowth);
end
end

% ============================================================ scope
function sc = i_prepareScope(model, GSS, relocate, comps, defaultCompartment, minGrowth, transportable, verbose) %#ok<INUSD>
biomassIdx = find(model.c ~= 0, 1);
sc.biomassId = model.rxns{biomassIdx};

if isempty(minGrowth)
    sol = solveLP(model);
    if isempty(sol.f) || abs(sol.f) <= 0
        error('RAVEN:badInput','the draft model does not grow; pass minGrowth.');
    end
    minGrowth = 0.1 * abs(sol.f);
end
sc.minGrowth = minGrowth;

% each reaction's single compartment ('' if boundary/multi-compartment)
sc.rxnComp = i_reactionCompartments(model);

toRelocate = ismember(model.rxns, relocate);
isBoundary = (sum(model.S ~= 0, 1)' == 1);
movable = toRelocate & ~isBoundary & ~cellfun(@isempty, sc.rxnComp);
movable(biomassIdx) = false;
sc.movIdx = find(movable);
sc.pinIdx = find(~movable);

% genes on movable reactions that are scored
[gInGSS, gssRow] = ismember(model.genes, GSS.genes);
geneOnMov = any(model.rxnGeneMat(sc.movIdx, :) ~= 0, 1)';
scope = gInGSS & geneOnMov;
sc.geneIdx = find(scope);
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

% base metabolite = metabolite name; transportable base names
sc.metNames = model.metNames;
movMetMask = any(model.S(:, sc.movIdx) ~= 0, 2);
if isnumeric(transportable) && isempty(transportable)
    sc.isTransp = movMetMask;                                  % default: all movable bases
else
    sc.isTransp = movMetMask & ismember(model.metNames, transportable);
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

pR=[];pC=[];pV=[];
for k=1:nMov
    for ci=1:nC; pR(end+1)=k; pC(end+1)=xCol(k,ci); pV(end+1)=1; end %#ok<AGROW>
end
A_place=sparse(pR,pC,pV,nMov,nVar); b_place=ones(nMov,1);

cR=[];cC=[];cV=[];row=0;
for k=1:nMov
    gOn=find(model.rxnGeneMat(movIdx(k),:)~=0);
    for g=gOn
        gi=find(geneIdx==g,1); if isempty(gi); continue; end
        for ci=1:nC
            row=row+1;
            cR(end+1)=row; cC(end+1)=xCol(k,ci); cV(end+1)=1;   %#ok<AGROW>
            cR(end+1)=row; cC(end+1)=yCol(gi,ci); cV(end+1)=-1; %#ok<AGROW>
        end
    end
end
A_couple=sparse(cR,cC,cV,row,nVar); b_couple=zeros(row,1);

hR=[];hC=[];hV=[];row=0;
for gi=1:nGene
    g=geneIdx(gi);
    rOfG=find(model.rxnGeneMat(movIdx,g)~=0)';
    for ci=1:nC
        row=row+1; hR(end+1)=row; hC(end+1)=yCol(gi,ci); hV(end+1)=1; %#ok<AGROW>
        for k=rOfG; hR(end+1)=row; hC(end+1)=xCol(k,ci); hV(end+1)=-1; end %#ok<AGROW>
    end
end
A_has=sparse(hR,hC,hV,row,nVar); b_has=zeros(row,1);

g1R=[];g1C=[];g1V=[];
for gi=1:nGene
    for ci=1:nC; g1R(end+1)=gi; g1C(end+1)=yCol(gi,ci); g1V(end+1)=1; end %#ok<AGROW>
end
A_gene1=sparse(g1R,g1C,g1V,nGene,nVar); b_gene1=ones(nGene,1);

fR=[];fC=[];nForce=0; fkeys=keys(forced);
for i=1:numel(fkeys)
    k=fkeys{i};
    nForce=nForce+1; fR(end+1)=nForce; fC(end+1)=xCol(k,forced(k)); %#ok<AGROW>
end
A_force=sparse(fR,fC,ones(1,nForce),nForce,nVar); b_force=ones(nForce,1);

coR=[];coC=[];coV=[];row=0;
for gi=1:numel(groups)
    mem=groups{gi};
    for j=1:numel(mem)-1
        a=mem(j); b=mem(j+1);
        for ci=1:nC
            row=row+1;
            coR(end+1)=row; coC(end+1)=xCol(a,ci); coV(end+1)=1;  %#ok<AGROW>
            coR(end+1)=row; coC(end+1)=xCol(b,ci); coV(end+1)=-1; %#ok<AGROW>
        end
    end
end
A_colo=sparse(coR,coC,coV,row,nVar); b_colo=zeros(row,1);

c = zeros(nVar,1);
for gi=1:nGene
    for ci=1:nC; c(yCol(gi,ci)) = score(gi,ci) - multiPen; end
end

prob.A = [A_place; A_couple; A_has; A_gene1; A_force; A_colo];
prob.a = prob.A;
prob.b = [b_place; b_couple; b_has; b_gene1; b_force; b_colo];
prob.csense = [repmat('E',1,nMov), repmat('L',1,numel(b_couple)), ...
               repmat('L',1,numel(b_has)), repmat('G',1,nGene), ...
               repmat('E',1,nForce), repmat('E',1,numel(b_colo))];
prob.c = -c; prob.osense = 1;
prob.lb = zeros(nVar,1); prob.ub = ones(nVar,1);
prob.vartype = repmat('B',1,nVar);

params.intTol = 1e-9; params.TimeLimit = 1000;
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

nMet = numel(model.mets);
usedComp = cell(nMet,1); usedRxn = cell(nMet,1);
allRxn = [sc.movIdx; sc.pinIdx];
for r = allRxn'
    comp = placedComp(r);
    if comp == 0; continue; end     % multi-compartment reaction: bridges pools
    for mm = find(model.S(:,r)~=0)'
        if sc.isTransp(mm); continue; end
        usedComp{mm}(end+1) = comp; usedRxn{mm}(end+1) = r;
    end
end

for mm=1:nMet
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
function [trMet, trComp] = i_splitTransports(model, sc, comps, defaultCompartment, placeIdx, relaxed)
defC = find(strcmp(comps, defaultCompartment), 1);
placedComp = i_placedCompartments(model, sc, comps, placeIdx);
allRxn = [sc.movIdx; sc.pinIdx];
nMet = numel(model.mets);
used = cell(nMet,1);
for r = allRxn'
    comp = placedComp(r);
    if comp == 0; continue; end
    for mm = find(model.S(:,r)~=0)'
        if sc.isTransp(mm) || (isKey(relaxed,mm))
            used{mm}(end+1) = comp;
        end
    end
end
trMet=[]; trComp=[];
for mm=1:nMet
    cs = unique(used{mm});
    if numel(cs) <= 1; continue; end
    for ci = sort(cs)
        if ci ~= defC; trMet(end+1)=mm; trComp(end+1)=ci; end %#ok<AGROW>
    end
end
trMet=trMet(:); trComp=trComp(:);
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
function [trMet, trComp] = i_usableTransports(outModel, trMet, trComp, comps, biomassId, growthConditions)
% Keep only transports that can carry flux in at least one certification
% medium (sound: an unusable reaction's removal cannot change any of those
% FBAs).
trId = cell(numel(trMet),1);
for i=1:numel(trMet)
    trId{i} = ['tr_' num2str(i-1) '_' comps{trComp(i)}];
end
media = {[]};
if ~isempty(growthConditions) && isstruct(growthConditions)
    for i=1:numel(growthConditions); media{end+1} = growthConditions(i).medium; end %#ok<AGROW>
end
keep = false(numel(trMet),1);
for mi=1:numel(media)
    m = outModel;
    if ~isempty(media{mi}); m = i_applyMedium(m, media{mi}); end
    bIdx = find(strcmp(m.rxns, biomassId), 1);
    m.c = zeros(numel(m.rxns),1); m.c(bIdx) = 1;
    idx = zeros(numel(trId),1);
    for i=1:numel(trId); j=find(strcmp(m.rxns,trId{i}),1); if ~isempty(j); idx(i)=j; end; end
    valid = idx>0;
    fl = false(numel(trId),1);
    fl(valid) = haveFlux(m, 'rxns', idx(valid));
    keep = keep | fl;
end
trMet = trMet(keep); trComp = trComp(keep);
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
function outModel = i_applyAssignment(model, sc, comps, defaultCompartment, placeIdx, trMet, trComp)
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
for k = 1:numel(trMet)
    outModel = i_addTransport(outModel, model, trMet(k), comps{trComp(k)}, defaultCompartment, k-1);
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
