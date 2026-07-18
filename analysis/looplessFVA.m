function [minFlux, maxFlux] = looplessFVA(model, rxns, minGrowth)
% looplessFVA  Loop-free flux variability analysis.
%
% Minimum and maximum flux of each listed reaction over the steady-state flux
% space, excluding any flux that could exist only around a thermodynamically
% infeasible internal loop. Uses the loop-law formulation of Schellenberger et
% al. (2011), Biophys J 100:544: a metabolite-potential vector mu assigns each
% internal reaction an energy G = S'*mu that must oppose its flux direction, so
% no internal cycle can carry net flux. Each reaction's min/max is then a MILP.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
% rxns : cell or double
%     reaction ids or indexes to analyse (default all reactions).
% minGrowth : double
%     if given, the model objective is held at or above this value while the
%     variability is computed (default []: the objective is not constrained).
%
% Returns
% -------
% minFlux : double
%     the loop-free minimum flux of each reaction in rxns.
% maxFlux : double
%     the loop-free maximum flux of each reaction in rxns.
%
% See also
% --------
% cycleFreeFlux, getAllowedBounds, solveLP

if nargin < 2 || isempty(rxns)
    idx = (1:numel(model.rxns))';
else
    idx = getIndexes(model, rxns, 'rxns');
end
if nargin < 3; minGrowth = []; end

M = 1000;
nRxns = numel(model.rxns);
nMets = numel(model.mets);
% Internal reactions get the loop law; exchange reactions cannot be in an
% internal loop, so they are exempt.
[~, exchIdx] = getExchangeRxns(model);
isExch = false(nRxns,1); isExch(exchIdx) = true;
intIdx = find(~isExch);
nInt = numel(intIdx);

% Variables: v (nRxns), mu (nMets), z (nInt binary).
oV = 0; oMu = nRxns; oZ = nRxns + nMets; nVar = nRxns + nMets + nInt;

% --- constraints ---
% S v = 0
A_S = [model.S, sparse(nMets, nMets + nInt)];
% per-internal loop-law rows
rA = []; cA = []; vA = []; bDir = []; senseDir = ''; row = 0;
for k = 1:nInt
    gi = intIdx(k); zc = oZ + k;
    lbi = model.lb(gi); ubi = model.ub(gi);
    % (a) v_gi - ub*z <= 0
    row=row+1; rA(end+1)=row; cA(end+1)=gi; vA(end+1)=1; rA(end+1)=row; cA(end+1)=zc; vA(end+1)=-ubi; bDir(end+1)=0; senseDir(end+1)='L'; %#ok<AGROW>
    % (b) v_gi + lb*z >= lb
    row=row+1; rA(end+1)=row; cA(end+1)=gi; vA(end+1)=1; rA(end+1)=row; cA(end+1)=zc; vA(end+1)=lbi; bDir(end+1)=lbi; senseDir(end+1)='G'; %#ok<AGROW>
    % (c) S(:,gi)'mu + M*z <= M-1
    mrows = find(model.S(:,gi));
    for mm=mrows'; row0=row+1; rA(end+1)=row0; cA(end+1)=oMu+mm; vA(end+1)=model.S(mm,gi); end %#ok<AGROW>
    row=row+1; rA(end+1)=row; cA(end+1)=zc; vA(end+1)=M; bDir(end+1)=M-1; senseDir(end+1)='L'; %#ok<AGROW>
    % (d) S(:,gi)'mu + M*z >= 1
    for mm=mrows'; row0=row+1; rA(end+1)=row0; cA(end+1)=oMu+mm; vA(end+1)=model.S(mm,gi); end %#ok<AGROW>
    row=row+1; rA(end+1)=row; cA(end+1)=zc; vA(end+1)=M; bDir(end+1)=1; senseDir(end+1)='G'; %#ok<AGROW>
end
A_dir = sparse(rA, cA, vA, row, nVar);

A = [A_S; A_dir];
b = [zeros(nMets,1); bDir(:)];
csense = [repmat('E',1,nMets), senseDir];

% objective (biomass) floor
if ~isempty(minGrowth)
    objIdx = find(model.c ~= 0);
    for o=objIdx'
        A(end+1,:) = sparse(1, oV+o, 1, 1, nVar); %#ok<AGROW>
        b(end+1,1) = minGrowth; csense(end+1) = 'G'; %#ok<AGROW>
    end
end

lb = [model.lb(:); -M*ones(nMets,1); zeros(nInt,1)];
ub = [model.ub(:);  M*ones(nMets,1);  ones(nInt,1)];
vartype = [repmat('C',1,nRxns+nMets), repmat('B',1,nInt)];

prob.a = A; prob.A = A; prob.b = b; prob.csense = csense;
prob.lb = lb; prob.ub = ub; prob.vartype = vartype; prob.osense = 1;
params.intTol = 1e-9; params.TimeLimit = 1000;

minFlux = zeros(numel(idx),1); maxFlux = zeros(numel(idx),1);
for k = 1:numel(idx)
    r = idx(k);
    c = zeros(nVar,1); c(oV+r) = 1;
    prob.c = -c;                       % maximise v_r
    sol = optimizeProb(prob, params, false);
    if checkSolution(sol); maxFlux(k) = sol.full(oV+r); end
    prob.c = c;                        % minimise v_r
    sol = optimizeProb(prob, params, false);
    if checkSolution(sol); minFlux(k) = sol.full(oV+r); end
end
end
