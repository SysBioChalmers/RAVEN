function [outModel, placement, addedTransports, exitFlag, report] = assignCompartments(model, GSS, reactionsToRelocate, varargin)
% assignCompartments  Assign reactions to compartments, certified by growth.
%
% Deterministic alternative to predictLocalization. Places the requested
% reactions into the compartments named by localization scores (GSS) and then
% verifies, by an actual FBA on the built model, that the objective (biomass)
% is still producible. The work is done in three phases:
%
%   1. Place   - a flux-free MILP that maximises the localization score.
%                Its only variables are placement binaries (reaction ->
%                compartment, gene -> compartment); it has no flux variables
%                and no growth constraint, so a placement can never harvest a
%                compartment's score through leaked flux.
%   2. Repair  - a solver-free fixpoint that keeps the placement connected: a
%                non-transportable metabolite split across compartments forces
%                its reactions to co-locate; transportable splits get passive
%                transports through the default compartment.
%   3. Certify - the placement is materialised into a compartmentalised model
%                and that exact model is solved. exitFlag reflects whether it
%                reached the growth floor, so the certificate is the model
%                that is returned, never a placement the model cannot grow.
%
% Mono-localization: each reaction is placed in exactly one compartment; a
% gene still spans compartments when it catalyses reactions placed in
% different ones.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model with an objective (model.c) and grRules. Multiple
%     compartments are merged before placement.
% GSS : struct
%     gene scoring structure (genes, compartments, scores) as from parseScores.
%     GSS.compartments are the target compartment labels and must include
%     defaultCompartment.
% reactionsToRelocate : cell
%     reaction ids to (re)place. Boundary reactions and the objective reaction
%     are always pinned.
%
% Name-Value Arguments
% --------------------
% defaultCompartment : char
%     compartment that transports route through (usually cytosol); must be in
%     GSS.compartments.
% multiCompartmentPenalty : double (default 0.5)
%     score cost per extra compartment a gene ends up in.
% minGrowth : double (default [] = 10%% of the unconstrained optimum)
%     required objective flux for certification.
% transportable : cell (default [] = all)
%     metabolite ids (merged-model names) that may receive transports.
%     Restricting it forces functionality-driven placement of reactions whose
%     substrates are confined.
% growthConditions : struct array (default [])
%     extra media the placement must also grow on, each with fields: name,
%     medium (struct of exchangeRxnId -> max uptake), and minGrowth.
% maxRounds : double (default 8)
%     maximum place/repair rounds before returning the best placement found.
% transportCost : double (default 0.5)
%     accepted for backward compatibility but not used: placement is driven by
%     the score alone, and transports are a structural consequence of repair.
% verbose : logical (default true)
%
% Returns
% -------
% outModel : struct
%     the compartmentalised model.
% placement : struct
%     .rxns and .compartment (the assigned compartment id per relocated
%     reaction).
% addedTransports : struct
%     .mets and .compartment for the transports repair added.
% exitFlag : double
%     1 = certified (the built model reaches the growth floor on every
%     medium), -1 = could not place or the built model did not certify.
% report : struct
%     .certified (logical), .status (char), and .growths (struct of
%     medium -> objective flux), so a failed certification reports the real
%     growth rather than hiding it.
%
% See also
% --------
% predictLocalization, parseScores, gapFillMILP

p = parseRAVENargs(varargin, {'defaultCompartment',[]; 'transportCost',0.5; ...
    'multiCompartmentPenalty',0.5; 'minGrowth',[]; 'transportable',[]; ...
    'growthConditions',[]; 'maxRounds',8; 'verbose',true});
defaultCompartment = char(p.defaultCompartment);
multiPen = p.multiCompartmentPenalty;
minGrowth = p.minGrowth;
growthConditions = p.growthConditions;
maxRounds = p.maxRounds;
verbose = p.verbose;

outModel = model; placement = struct('rxns',{{}},'compartment',{{}});
addedTransports = struct('mets',{{}},'compartment',{{}}); exitFlag = -1;
report = struct('certified',false,'status','not_solved','growths',struct());

if all(model.c == 0)
    error('RAVEN:badInput','model has no objective (set model.c).');
end
comps = GSS.compartments(:);
if ~ismember(defaultCompartment, comps)
    error('RAVEN:badInput','defaultCompartment ''%s'' not in GSS.compartments.', defaultCompartment);
end

% Merge to a single compartment so every metabolite has one identity;
% placement and transports are then expressed relative to defaultCompartment.
% After merging, treat that single compartment AS the default: pinned
% reactions are considered to sit there, and placed reactions and transports
% are expressed relative to it.
if numel(model.comps) > 1
    model = mergeCompartments(model, true, true);
    model.comps = {defaultCompartment};
    if isfield(model,'compNames'); model.compNames = {defaultCompartment}; end
    if isfield(model,'compOutside'); model.compOutside = {''}; end
    model.metComps = ones(numel(model.mets),1);
end
biomassIdx = find(model.c ~= 0, 1);

if isempty(minGrowth)
    sol = solveLP(model);
    if isempty(sol.f) || sol.f <= 0
        error('RAVEN:badInput','merged model does not grow; pass minGrowth.');
    end
    minGrowth = 0.1 * abs(sol.f);
end

% ---- scope: movable vs pinned ----
isBoundary = (sum(model.S ~= 0, 1)' == 1);   % one-metabolite reactions
relocate = ismember(model.rxns, reactionsToRelocate);
movableMask = relocate & ~isBoundary;
movableMask(biomassIdx) = false;
movIdx = find(movableMask);
pinIdx = find(~movableMask);
nMov = numel(movIdx); nC = numel(comps);
defC = find(strcmp(comps, defaultCompartment));

% Genes in scope: those on movable reactions that have a row in GSS
[gInGSS, gssRow] = ismember(model.genes, GSS.genes);
geneOnMov = any(model.rxnGeneMat(movIdx, :) ~= 0, 1)';
geneIdx = find(gInGSS & geneOnMov);
nGene = numel(geneIdx);
score = zeros(nGene, nC);
for gi = 1:nGene
    score(gi, :) = GSS.scores(gssRow(geneIdx(gi)), :);
end

% Transportable metabolites (a base metabolite may be moved between
% compartments by a passive transport during repair).
movMetMask = any(model.S(:, movIdx) ~= 0, 2);
if isnumeric(p.transportable) && isempty(p.transportable)
    transpMet = find(movMetMask);                                   % default: all
else
    transpMet = find(movMetMask & ismember(model.mets, p.transportable));
end
isTransp = false(numel(model.mets),1); isTransp(transpMet) = true;

if verbose
    fprintf('assignCompartments: %d movable, %d pinned reactions, %d genes, %d compartments.\n', ...
        nMov, numel(pinIdx), nGene, nC);
end

% ---- place / repair rounds ----
% forced(mi) = compartment index a movable reaction is pinned to (0 = free);
% groups is a cell array of movable-index vectors that must share a
% compartment. Both grow monotonically until the confinement fixpoint holds.
forced = zeros(nMov,1);
groups = {};
seen = {};
placeIdx = [];
for round = 1:maxRounds
    placeIdx = i_placementMaster(model, movIdx, geneIdx, score, multiPen, ...
        nC, forced, groups, verbose);
    if isempty(placeIdx)
        report.status = 'infeasible';
        if verbose; fprintf('assignCompartments: placement MILP infeasible.\n'); end
        return;
    end
    [newForced, newGroups] = i_diagnoseConfinement(model, movIdx, pinIdx, ...
        placeIdx, defC, isTransp);
    % Apply only genuinely new tightenings; stop at the fixpoint.
    changed = false;
    for k = 1:nMov
        if newForced(k) ~= 0 && forced(k) == 0
            forced(k) = newForced(k); changed = true;
        end
    end
    for gi = 1:numel(newGroups)
        if ~i_hasGroup(groups, newGroups{gi})
            groups{end+1} = newGroups{gi}; changed = true; %#ok<AGROW>
        end
    end
    sig = i_signature(forced, groups);
    if ~changed || any(strcmp(seen, sig))
        break;
    end
    seen{end+1} = sig; %#ok<AGROW>
end

% ---- placement + transports ----
placeRxns = model.rxns(movIdx);
placeComp = comps(placeIdx);
placement.rxns = placeRxns(:); placement.compartment = placeComp(:);

[trMet, trComp] = i_splitTransports(model, movIdx, pinIdx, placeIdx, defC, isTransp, nC);
addedTransports.mets = model.mets(trMet);
addedTransports.compartment = comps(trComp);

% ---- materialise + certify ----
outModel = i_applyAssignment(model, placement, addedTransports, comps, defaultCompartment);
[certified, growths] = i_certify(outModel, biomassIdx, model.rxns{biomassIdx}, minGrowth, growthConditions);
report.certified = certified;
report.growths = growths;
if certified
    exitFlag = 1; report.status = 'certified';
else
    exitFlag = -1; report.status = 'uncertified';
    if verbose
        fprintf('assignCompartments: placement did not certify (growth %.4g < %.4g).\n', ...
            growths.primary, minGrowth);
    end
end
end

% ------------------------------------------------------------- placement MILP
function placeIdx = i_placementMaster(model, movIdx, geneIdx, score, multiPen, nC, forced, groups, verbose)
% Flux-free score MILP. Variables: x[mov,c] and y[gene,c], both binary.
% Returns the chosen compartment index per movable reaction, or [] if
% infeasible.
nMov = numel(movIdx); nGene = numel(geneIdx);
oX = 0;      nX = nMov*nC;
oY = oX+nX;  nY = nGene*nC;
nVar = oY+nY;
xCol = @(mi,ci) oX + (mi-1)*nC + ci;
yCol = @(gi,ci) oY + (gi-1)*nC + ci;

% placement: sum_c x[r,c] = 1
pR=[];pC=[];pV=[];
for k=1:nMov
    for ci=1:nC; pR(end+1)=k; pC(end+1)=xCol(k,ci); pV(end+1)=1; end %#ok<AGROW>
end
A_place=sparse(pR,pC,pV,nMov,nVar); b_place=ones(nMov,1);

% coupling: x[r,c] - y[g,c] <= 0 for each gene g on reaction r
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

% has: y[g,c] - sum_{r of g} x[r,c] <= 0
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

% gene1: sum_c y[g,c] >= 1
g1R=[];g1C=[];g1V=[];
for gi=1:nGene
    for ci=1:nC; g1R(end+1)=gi; g1C(end+1)=yCol(gi,ci); g1V(end+1)=1; end %#ok<AGROW>
end
A_gene1=sparse(g1R,g1C,g1V,nGene,nVar); b_gene1=ones(nGene,1);

% forced pins: x[r,c_force] = 1
fR=[];fC=[];nForce=0;
for k=1:nMov
    if forced(k)~=0
        nForce=nForce+1; fR(end+1)=nForce; fC(end+1)=xCol(k,forced(k)); %#ok<AGROW>
    end
end
A_force=sparse(fR,fC,ones(1,nForce),nForce,nVar); b_force=ones(nForce,1);

% co-location groups: x[a,c] - x[b,c] = 0 for consecutive members
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

% objective: max sum score*y - multiPen*sum y, with a small score-consistent
% tie-break on the reaction placement. The gene reward on y decides which
% compartment each gene occupies, but leaves the reaction binaries x
% underdetermined (any placement consistent with the gene assignment is
% optimal); without the tie-break the solver picks one arbitrarily. The tiny
% reward on x pulls each reaction into the compartment its own genes score
% highest, without being large enough to change the gene assignment.
xTie = 1e-3;
c = zeros(nVar,1);
for gi=1:nGene
    for ci=1:nC; c(yCol(gi,ci)) = score(gi,ci) - multiPen; end
end
for k=1:nMov
    gOn = find(model.rxnGeneMat(movIdx(k),:)~=0);
    for g=gOn
        gi=find(geneIdx==g,1); if isempty(gi); continue; end
        for ci=1:nC; c(xCol(k,ci)) = c(xCol(k,ci)) + xTie*score(gi,ci); end
    end
end

prob.A = [A_place; A_couple; A_has; A_gene1; A_force; A_colo];
prob.a = prob.A;
prob.b = [b_place; b_couple; b_has; b_gene1; b_force; b_colo];
prob.csense = [repmat('E',1,nMov), repmat('L',1,numel(b_couple)), ...
               repmat('L',1,numel(b_has)), repmat('G',1,nGene), ...
               repmat('E',1,nForce), repmat('E',1,numel(b_colo))];
prob.c = -c;            % optimizeProb minimises; we maximise c
prob.osense = 1;
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

% ------------------------------------------------------ confinement diagnosis
function [forced, groups] = i_diagnoseConfinement(model, movIdx, pinIdx, placeIdx, defC, isTransp)
% Solver-free: find non-transportable metabolites split across compartments
% and decide, per split, whether to force the movers to the pinned
% compartment or to co-locate them. Pinned reactions all live in the default
% compartment (the model was merged first).
nMov = numel(movIdx);
forced = zeros(nMov,1);
groups = {};

% used{metRow} maps compartment index -> list of movable local indices.
% pinnedHere(metRow) is true if a pinned reaction touches the metabolite.
nMet = numel(model.mets);
usedComp = cell(nMet,1);
usedMov  = cell(nMet,1);
pinnedHere = false(nMet,1);
for k=1:nMov
    mrows = find(model.S(:,movIdx(k))~=0)';
    ci = placeIdx(k);
    for mm=mrows
        if isTransp(mm); continue; end     % transportable: handled by a transport
        usedComp{mm}(end+1) = ci;
        usedMov{mm}(end+1) = k;
    end
end
for k=1:numel(pinIdx)
    mrows = find(model.S(:,pinIdx(k))~=0)';
    for mm=mrows
        if isTransp(mm); continue; end
        usedComp{mm}(end+1) = defC;
        pinnedHere(mm) = true;
    end
end

for mm=1:nMet
    cs = unique(usedComp{mm});
    if numel(cs) <= 1; continue; end       % lives in one compartment: fine
    movers = unique(usedMov{mm});
    if pinnedHere(mm)
        % A pinned reaction anchors this metabolite in the default
        % compartment; every movable toucher must join it there.
        for k=movers; forced(k) = defC; end
    else
        % Only movable reactions touch it: co-locate them and let the score
        % objective choose the shared compartment.
        groups{end+1} = movers(:)'; %#ok<AGROW>
    end
end
end

% ----------------------------------------------------------- split transports
function [trMet, trComp] = i_splitTransports(model, movIdx, pinIdx, placeIdx, defC, isTransp, nC)
% For each transportable base metabolite placed in more than one compartment,
% add a transport between the default compartment and each non-default one.
nMet = numel(model.mets);
compsUsed = cell(nMet,1);
for k=1:numel(movIdx)
    mrows = find(model.S(:,movIdx(k))~=0)';
    for mm=mrows; compsUsed{mm}(end+1) = placeIdx(k); end
end
for k=1:numel(pinIdx)
    mrows = find(model.S(:,pinIdx(k))~=0)';
    for mm=mrows; compsUsed{mm}(end+1) = defC; end
end
trMet=[]; trComp=[];
for mm=1:nMet
    if ~isTransp(mm); continue; end
    cs = unique(compsUsed{mm});
    if numel(cs) <= 1; continue; end
    for ci=cs
        if ci==defC; continue; end
        trMet(end+1)=mm; trComp(end+1)=ci; %#ok<AGROW>
    end
end
trMet=trMet(:); trComp=trComp(:);
end

% ------------------------------------------------------------- certification
function [certified, growths] = i_certify(outModel, biomassIdx, biomassId, minGrowth, growthConditions)
% Solve the materialised model and confirm it reaches minGrowth on the primary
% medium and on every extra growth condition.
tol = 1e-9;
growths = struct();
bIdx = find(strcmp(outModel.rxns, biomassId), 1);
if isempty(bIdx); bIdx = biomassIdx; end
sol = solveLP(outModel);
primary = 0;
if ~isempty(sol.f); primary = abs(sol.f); end
growths.primary = primary;
certified = primary >= minGrowth - tol;

if ~isempty(growthConditions) && isstruct(growthConditions)
    for i=1:numel(growthConditions)
        gc = growthConditions(i);
        m2 = i_applyMedium(outModel, gc.medium);
        s2 = solveLP(m2);
        g = 0; if ~isempty(s2.f); g = abs(s2.f); end
        growths.(matlab.lang.makeValidName(gc.name)) = g;
        certified = certified && (g >= gc.minGrowth - tol);
    end
end
end

function model = i_applyMedium(model, medium)
% Close all uptake, then open only the listed exchanges. medium is a struct
% of exchangeRxnId -> max uptake.
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

% ------------------------------------------------------------------ helpers
function tf = i_hasGroup(groups, g)
tf = false;
g = sort(g(:)');
for i=1:numel(groups)
    if isequal(sort(groups{i}(:)'), g); tf = true; return; end
end
end

function s = i_signature(forced, groups)
gp = '';
for i=1:numel(groups); gp = [gp '|' num2str(sort(groups{i}(:)'))]; end %#ok<AGROW>
s = [num2str(forced(:)') '#' gp];
end

% ------------------------------------------------------------------- apply
function outModel = i_applyAssignment(model, placement, addedTransports, comps, defaultCompartment)
% Build the compartmentalised model: relabel each placed reaction's metabolites
% into its compartment (creating per-compartment metabolite copies) and add the
% requested transports.
outModel = model;
for ci = 1:numel(comps)
    if ~ismember(comps{ci}, outModel.comps)
        outModel.comps{end+1,1} = comps{ci};
        if isfield(outModel,'compNames'); outModel.compNames{end+1,1} = comps{ci}; end
    end
end
for k = 1:numel(placement.rxns)
    r = find(strcmp(outModel.rxns, placement.rxns{k}), 1);
    comp = placement.compartment{k};
    srcR = find(strcmp(model.rxns,placement.rxns{k}),1);
    mrows = find(model.S(:, srcR) ~= 0); %#ok<*FNDSB>
    for mm = mrows'
        [outModel, newMet] = i_metInComp(outModel, model, mm, comp);
        outModel.S(newMet, r) = model.S(mm, srcR);
        if ~isequal(newMet, mm); outModel.S(mm, r) = 0; end
    end
end
for k = 1:numel(addedTransports.mets)
    outModel = i_addTransport(outModel, model, addedTransports.mets{k}, addedTransports.compartment{k}, defaultCompartment);
end
end

function [outModel, idx] = i_metInComp(outModel, model, srcMetRow, comp)
ci = find(strcmp(outModel.comps, comp), 1);
name = model.metNames{srcMetRow};
cand = find(strcmp(outModel.metNames, name) & outModel.metComps == ci, 1);
if ~isempty(cand); idx = cand; return; end
newId = [model.mets{srcMetRow} '_' comp];
outModel.mets{end+1,1} = newId;
outModel.metNames{end+1,1} = name;
outModel.metComps(end+1,1) = ci;
outModel.S(end+1, :) = 0;
if isfield(outModel,'b'); outModel.b(end+1,1) = 0; end
idx = numel(outModel.mets);
end

function outModel = i_addTransport(outModel, model, metId, comp, defaultCompartment)
srcRow = find(strcmp(model.mets, metId), 1);
if isempty(srcRow); return; end
[outModel, dRow] = i_metInComp(outModel, model, srcRow, defaultCompartment);
[outModel, cRow] = i_metInComp(outModel, model, srcRow, comp);
trId = ['tr_' model.metNames{srcRow} '_' comp];
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
