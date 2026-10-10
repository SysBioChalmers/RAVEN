function [outModel, placement, addedTransports, exitFlag] = assignCompartments(model, GSS, reactionsToRelocate, varargin)
% assignCompartments  Assign reactions to compartments by a functionality-constrained MILP.
%
% Deterministic alternative to predictLocalization: a single MILP places the requested reactions
% into compartments to agree with localization scores (GSS) WHILE keeping the model's objective
% (biomass) producible. Because functionality is a hard constraint, a reaction is placed against
% its own top score whenever the network needs it there (pathway coherence emerges from requiring
% flux, not from a heuristic). Mono-localization: each reaction is placed in exactly one
% compartment; a gene still spans compartments when it catalyses reactions placed in different ones.
%
% The model is merged to a single compartment first (like predictLocalization); the MILP then
% re-distributes the relocatable reactions across the compartments named in GSS, adding passive
% transports (defaultCompartment <-> c) where needed.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model with an objective (model.c) and grRules. Multiple compartments are merged.
% GSS : struct
%     gene scoring structure (genes, compartments, scores) as from parseScores. GSS.compartments
%     are the target compartment labels and must include defaultCompartment.
% reactionsToRelocate : cell
%     reaction ids to (re)place. Boundary reactions and the objective reaction are always pinned.
%
% Name-Value Arguments
% --------------------
% defaultCompartment : char
%     compartment that transports route through (usually cytosol); must be in GSS.compartments.
% transportCost : double (default 0.5)
%     cost per added inter-compartment transport.
% multiCompartmentPenalty : double (default 0.5)
%     cost per extra compartment a gene ends up in.
% minGrowth : double (default [] = 10%% of the unconstrained optimum)
%     required objective flux.
% transportable : cell (default [] = all)
%     metabolite ids (merged-model names) that may receive transports. Restricting it forces
%     functionality-driven placement of reactions whose substrates are confined.
% bigM : double (default 1000)
%     Big-M for flux gating.
% verbose : logical (default true)
%
% Returns
% -------
% outModel : struct
%     the compartmentalised model (objective still producible).
% placement : struct
%     .rxns and .compartment (the assigned compartment id per relocated reaction).
% addedTransports : struct
%     .mets and .compartment for the transports the MILP added.
% exitFlag : double
%     1 = optimal, -1 = infeasible/failed.
%
% See also
% --------
% predictLocalization, parseScores, gapFillMILP

p = parseRAVENargs(varargin, {'defaultCompartment',[]; 'transportCost',0.5; ...
    'multiCompartmentPenalty',0.5; 'minGrowth',[]; 'transportable',[]; 'bigM',1000; 'verbose',true});
defaultCompartment = char(p.defaultCompartment);
transportCost = p.transportCost;
multiPen = p.multiCompartmentPenalty;
minGrowth = p.minGrowth;
M = p.bigM;
verbose = p.verbose;

outModel = model; placement = struct('rxns',{{}},'compartment',{{}});
addedTransports = struct('mets',{{}},'compartment',{{}}); exitFlag = -1;

if all(model.c == 0)
    error('RAVEN:badInput','model has no objective (set model.c).');
end
comps = GSS.compartments(:);
if ~ismember(defaultCompartment, comps)
    error('RAVEN:badInput','defaultCompartment ''%s'' not in GSS.compartments.', defaultCompartment);
end

% Merge to a single compartment so every metabolite has one identity; relocations and transports
% are then expressed relative to defaultCompartment.
if numel(model.comps) > 1
    model = mergeCompartments(model, true, true);
end
nMet = numel(model.mets);
biomassIdx = find(model.c ~= 0, 1);

% Determine min growth from the merged (single-compartment) model
if isempty(minGrowth)
    sol = solveLP(model);
    if isempty(sol.f) || sol.f <= 0
        error('RAVEN:badInput','merged model does not grow; pass minGrowth.');
    end
    minGrowth = 0.1 * sol.f;
end

% ---- scope: movable vs pinned ----
isBoundary = (sum(model.S ~= 0, 1)' == 1);   % one-metabolite reactions
relocate = ismember(model.rxns, reactionsToRelocate);
movableMask = relocate & ~isBoundary;
movableMask(biomassIdx) = false;
movIdx = find(movableMask);
pinIdx = find(~movableMask);
nMov = numel(movIdx); nPin = numel(pinIdx); nC = numel(comps);
defC = find(strcmp(comps, defaultCompartment));

% Genes in scope: those on movable reactions that have a row in GSS
[gInGSS, gssRow] = ismember(model.genes, GSS.genes);
geneOnMov = any(model.rxnGeneMat(movIdx, :) ~= 0, 1)';
scopeMask = gInGSS & geneOnMov;
geneIdx = find(scopeMask);
nGene = numel(geneIdx);
score = zeros(nGene, nC);   % gene x compartment score, columns aligned to comps
for gi = 1:nGene
    score(gi, :) = GSS.scores(gssRow(geneIdx(gi)), :);
end

% Transportable metabolites (touched by movable reactions)
movMetMask = any(model.S(:, movIdx) ~= 0, 2);
if isnumeric(p.transportable) && isempty(p.transportable)
    transpMet = find(movMetMask);   % default ([]): all movable metabolites transportable
else
    transpMet = find(movMetMask & ismember(model.mets, p.transportable));  % a (possibly empty) set
end
% transport variables only for non-default compartments
trPairs = [];   % [metRow, compCol]
for ci = 1:nC
    if ci == defC; continue; end
    trPairs = [trPairs; [transpMet, repmat(ci, numel(transpMet), 1)]]; %#ok<AGROW>
end
nTr = size(trPairs, 1);

if verbose
    fprintf('assignCompartments: %d movable, %d pinned reactions, %d genes, %d compartments, %d transports.\n', ...
        nMov, nPin, nGene, nC, nTr);
end

% ---- variable layout (columns) ----
% fmove[mov,c]  : nMov*nC continuous
% fpin[pin]     : nPin    continuous
% ftr[tr]       : nTr     continuous
% x[mov,c]      : nMov*nC binary
% y[gene,c]     : nGene*nC binary
% t[tr]         : nTr     binary
oFmove = 0;            nFmove = nMov*nC;
oFpin  = oFmove+nFmove; nFpin = nPin;
oFtr   = oFpin+nFpin;   % nTr
oX     = oFtr+nTr;     nX = nMov*nC;
oY     = oX+nX;        nY = nGene*nC;
oT     = oY+nY;        % nTr
nVar   = oT+nTr;
fmoveCol = @(mi,ci) oFmove + (mi-1)*nC + ci;
xCol     = @(mi,ci) oX     + (mi-1)*nC + ci;
yCol     = @(gi,ci) oY     + (gi-1)*nC + ci;

% ---- bounds ----
lb = zeros(nVar,1); ub = zeros(nVar,1);
for k = 1:nMov
    r = movIdx(k);
    for ci = 1:nC
        lb(fmoveCol(k,ci)) = min(model.lb(r),0); ub(fmoveCol(k,ci)) = max(model.ub(r),0);
    end
end
lb(oFpin+(1:nPin)) = model.lb(pinIdx); ub(oFpin+(1:nPin)) = model.ub(pinIdx);
lb(oFtr+(1:nTr)) = -M; ub(oFtr+(1:nTr)) = M;
lb(oX+(1:nX)) = 0; ub(oX+(1:nX)) = 1;
lb(oY+(1:nY)) = 0; ub(oY+(1:nY)) = 1;
lb(oT+(1:nTr)) = 0; ub(oT+(1:nTr)) = 1;

% ---- node balance: S over (met, compartment) ----
% node row index = (metRow-1)*nC + ci
nNode = nMet*nC;
nodeRow = @(mr,ci) (mr-1)*nC + ci;
ri = []; ci_ = []; vv = [];
% movable: metabolites placed in c
for k = 1:nMov
    r = movIdx(k);
    mrows = find(model.S(:,r) ~= 0);
    for ci = 1:nC
        for mm = mrows'
            ri(end+1)=nodeRow(mm,ci); ci_(end+1)=fmoveCol(k,ci); vv(end+1)=model.S(mm,r); %#ok<AGROW>
        end
    end
end
% pinned: metabolites in defaultCompartment
for k = 1:nPin
    r = pinIdx(k);
    mrows = find(model.S(:,r) ~= 0);
    for mm = mrows'
        ri(end+1)=nodeRow(mm,defC); ci_(end+1)=oFpin+k; vv(end+1)=model.S(mm,r); %#ok<AGROW>
    end
end
% transports: -1 at (met,default), +1 at (met,c)
for k = 1:nTr
    mm = trPairs(k,1); cc = trPairs(k,2);
    ri(end+1)=nodeRow(mm,defC); ci_(end+1)=oFtr+k; vv(end+1)=-1; %#ok<AGROW>
    ri(end+1)=nodeRow(mm,cc);   ci_(end+1)=oFtr+k; vv(end+1)=+1; %#ok<AGROW>
end
A_node = sparse(ri, ci_, vv, nNode, nVar);
keepNode = any(A_node ~= 0, 2);          % drop empty nodes
A_node = A_node(keepNode, :);
nNodeKept = size(A_node,1);

% ---- flux gating: lb*x <= fmove <= ub*x ----
giR=[]; giC=[]; giV=[]; row=0; gUb=[]; gLb=[];
for k=1:nMov
    r=movIdx(k);
    for ci=1:nC
        row=row+1;  % ub: fmove - ub*x <= 0
        giR(end+1)=row; giC(end+1)=fmoveCol(k,ci); giV(end+1)=1; %#ok<AGROW>
        giR(end+1)=row; giC(end+1)=xCol(k,ci); giV(end+1)=-max(model.ub(r),0); %#ok<AGROW>
        gUb(end+1)=0; %#ok<AGROW>
    end
end
A_gub = sparse(giR,giC,giV,row,nVar); b_gub=zeros(row,1);
giR=[]; giC=[]; giV=[]; row=0;
for k=1:nMov
    r=movIdx(k);
    for ci=1:nC
        row=row+1;  % lb: fmove - lb*x >= 0
        giR(end+1)=row; giC(end+1)=fmoveCol(k,ci); giV(end+1)=1; %#ok<AGROW>
        giR(end+1)=row; giC(end+1)=xCol(k,ci); giV(end+1)=-min(model.lb(r),0); %#ok<AGROW>
    end
end
A_glb = sparse(giR,giC,giV,row,nVar); b_glb=zeros(row,1);

% ---- transport gating: -M*t <= ftr <= M*t ----
tR=[];tC=[];tV=[];
for k=1:nTr
    tR(end+1)=k; tC(end+1)=oFtr+k; tV(end+1)=1; tR(end+1)=k; tC(end+1)=oT+k; tV(end+1)=-M; %#ok<AGROW>
end
A_tub = sparse(tR,tC,tV,nTr,nVar);  % ftr - M*t <= 0
tR=[];tC=[];tV=[];
for k=1:nTr
    tR(end+1)=k; tC(end+1)=oFtr+k; tV(end+1)=1; tR(end+1)=k; tC(end+1)=oT+k; tV(end+1)=M; %#ok<AGROW>
end
A_tlb = sparse(tR,tC,tV,nTr,nVar);  % ftr + M*t >= 0

% ---- placement: sum_c x[mov,c] = 1 ----
pR=[];pC=[];pV=[];
for k=1:nMov
    for ci=1:nC; pR(end+1)=k; pC(end+1)=xCol(k,ci); pV(end+1)=1; end %#ok<AGROW>
end
A_place=sparse(pR,pC,pV,nMov,nVar); b_place=ones(nMov,1);

% ---- gene coupling x[r,c] <= y[g,c]; gene-has y<=sum x; gene1 sum_c y>=1 ----
cR=[];cC=[];cV=[];row=0;
for k=1:nMov
    gOn=find(model.rxnGeneMat(movIdx(k),:)~=0);
    for g=gOn
        gi=find(geneIdx==g,1); if isempty(gi); continue; end
        for ci=1:nC
            row=row+1; cR(end+1)=row; cC(end+1)=xCol(k,ci); cV(end+1)=1;
            cR(end+1)=row; cC(end+1)=yCol(gi,ci); cV(end+1)=-1; %#ok<AGROW>
        end
    end
end
A_couple=sparse(cR,cC,cV,row,nVar); b_couple=zeros(row,1);
% gene-has: y[g,c] - sum_{r of g} x[r,c] <= 0
hR=[];hC=[];hV=[];row=0;
for gi=1:nGene
    g=geneIdx(gi);
    rOfG=find(model.rxnGeneMat(movIdx,g)~=0)';  % indices into movIdx
    for ci=1:nC
        row=row+1; hR(end+1)=row; hC(end+1)=yCol(gi,ci); hV(end+1)=1;
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

% ---- growth: fpin[biomass] >= minGrowth ----
pinBio = find(pinIdx==biomassIdx,1);
A_grow = sparse(1, oFpin+pinBio, 1, 1, nVar); b_grow = minGrowth;

% ---- objective: max sum score*y - multiPen*sum y - transportCost*sum t ----
c = zeros(nVar,1);
for gi=1:nGene
    for ci=1:nC; c(yCol(gi,ci)) = c(yCol(gi,ci)) + score(gi,ci) - multiPen; end
end
if isscalar(transportCost); tcost = repmat(transportCost,nTr,1); else; tcost = transportCost(trPairs(:,1)); end
c(oT+(1:nTr)) = c(oT+(1:nTr)) - tcost(:);

% ---- assemble ----
prob.A = [A_node; A_gub; A_glb; A_tub; A_tlb; A_place; A_couple; A_has; A_gene1; A_grow];
prob.a = prob.A;
prob.b = [zeros(nNodeKept,1); b_gub; b_glb; zeros(nTr,1); zeros(nTr,1); b_place; b_couple; b_has; b_gene1; b_grow];
prob.csense = [repmat('E',1,nNodeKept), repmat('L',1,numel(b_gub)), repmat('G',1,numel(b_glb)), ...
               repmat('L',1,nTr), repmat('G',1,nTr), repmat('E',1,nMov), ...
               repmat('L',1,numel(b_couple)), repmat('L',1,numel(b_has)), repmat('G',1,nGene), 'G'];
prob.c = -c;          % optimizeProb minimises; we maximise c
prob.osense = 1;
prob.lb = lb; prob.ub = ub;
prob.vartype = repmat('C',1,nVar);
prob.vartype(oX+(1:nX)) = 'B';
prob.vartype(oY+(1:nY)) = 'B';
prob.vartype(oT+(1:nTr)) = 'B';

params.intTol = 1e-9; params.TimeLimit = 1000;
sol = optimizeProb(prob, params, verbose);
if ~checkSolution(sol)
    if verbose; fprintf('assignCompartments: MILP infeasible or failed.\n'); end
    return;
end

% ---- extract ----
xval = sol.full(oX+(1:nX));
placeRxns = {}; placeComp = {};
for k=1:nMov
    [~,ci] = max(xval((k-1)*nC + (1:nC)));
    placeRxns{end+1} = model.rxns{movIdx(k)}; placeComp{end+1} = comps{ci}; %#ok<AGROW>
end
placement.rxns = placeRxns(:); placement.compartment = placeComp(:);
tval = sol.full(oT+(1:nTr));
sel = find(tval > 0.5);
addedTransports.mets = model.mets(trPairs(sel,1));
addedTransports.compartment = comps(trPairs(sel,2));
exitFlag = 1;

outModel = i_applyAssignment(model, placement, addedTransports, comps, defaultCompartment);
end

% ------------------------------------------------------------------- apply
function outModel = i_applyAssignment(model, placement, addedTransports, comps, defaultCompartment)
% Build the compartmentalised model: relabel each placed reaction's metabolites into its
% compartment (creating per-compartment metabolite copies) and add the requested transports.
outModel = model;
% ensure all target compartments exist
for ci = 1:numel(comps)
    if ~ismember(comps{ci}, outModel.comps)
        outModel.comps{end+1,1} = comps{ci};
        if isfield(outModel,'compNames'); outModel.compNames{end+1,1} = comps{ci}; end
    end
end
% metabolite (baseName, compartment) -> index, seeded with existing
metKey = strcat(model.metNames, '###', model.comps(model.metComps));
for k = 1:numel(placement.rxns)
    r = find(strcmp(outModel.rxns, placement.rxns{k}), 1);
    comp = placement.compartment{k};
    mrows = find(model.S(:, find(strcmp(model.rxns,placement.rxns{k}),1)) ~= 0); %#ok<*FNDSB>
    for mm = mrows'
        [outModel, newMet] = i_metInComp(outModel, model, mm, comp);
        outModel.S(newMet, r) = model.S(mm, find(strcmp(model.rxns,placement.rxns{k}),1));
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
if ~isempty(cand); idx = cand; return; end   % reuse existing per-compartment metabolite
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
