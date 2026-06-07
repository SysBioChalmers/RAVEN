function [prob, condInfo] = splitProbForConditioning(prob, maxRatio)
% splitProbForConditioning
%   Reformulates an LP problem to improve its numerical conditioning by
%   splitting columns whose coefficients span more than a given ratio. In a
%   metabolic model a column corresponds to a reaction, and a wide spread of
%   stoichiometric coefficients within one reaction (e.g. a biomass reaction
%   with both major substrates and trace cofactors) cannot be removed by the
%   diagonal row/column scaling that LP solvers apply, and may lead to
%   unreliable solutions or spurious infeasibility.
%
%   Each offending small coefficient is replaced by a chain consisting of an
%   auxiliary "star" metabolite (a new equality constraint row) and one or
%   more conversion variables (new columns), such that no individual column
%   has a coefficient ratio exceeding maxRatio. The feasible region of the
%   original problem is preserved exactly: the new rows are steady-state (=0)
%   constraints and the new variables are fully determined by them, so the
%   values of all original variables are unchanged.
%
%   The function is meant to be called inside optimizeProb, on the COBRA-style
%   problem struct, immediately before the problem is passed to the solver.
%   The auxiliary rows are appended at the bottom of prob.A and the auxiliary
%   columns at the right, so the indices of all original constraints and
%   variables are preserved. This allows the augmentation to be stripped from
%   the solver output (see condInfo).
%
% Input:
%   prob        COBRA-style LP problem struct, with at least the fields A, b,
%               c, lb, ub, csense and vartype
%   maxRatio    maximum allowed ratio between the largest and smallest
%               absolute coefficient within any single column
%
% Output:
%   prob        the (possibly) augmented problem struct. Identical to the
%               input if no column needed splitting
%   condInfo    struct describing the augmentation:
%       applied     true if any column was split
%       nSplits     number of auxiliary metabolite/conversion pairs added
%       nOrigRows   number of constraint rows before augmentation
%       nOrigCols   number of variables (columns) before augmentation
%
% Usage: [prob, condInfo] = splitProbForConditioning(prob, maxRatio)

[nOrigRows, nOrigCols] = size(prob.A);
condInfo.applied   = false;
condInfo.nSplits   = 0;
condInfo.nOrigRows = nOrigRows;
condInfo.nOrigCols = nOrigCols;

%Detection uses a small slack so that reactions that are only marginally over
%the threshold are left alone. This also makes the transformation idempotent:
%a freshly split column sits exactly at maxRatio and is not split again.
detectRatio = maxRatio * 1.01;

%Work on triplets to avoid repeated resizing of the sparse matrix
[I, J, V] = find(prob.A);
absV = abs(V);

%Largest/smallest absolute coefficient and number of nonzeros per column
cMax   = accumarray(J, absV, [nOrigCols 1], @max, 0);
cMin   = accumarray(J, absV, [nOrigCols 1], @min, inf);
nnzCol = accumarray(J, 1,    [nOrigCols 1], @sum, 0);

%Single-coefficient columns (e.g. exchange reactions or slack variables) have
%an undefined ratio and are never split
badCols = find(nnzCol >= 2 & (cMax ./ cMin) > detectRatio);
if isempty(badCols)
    return
end

%Flag the entries that are too small within their (bad) column. These are
%removed from the original matrix and re-introduced via the auxiliary chain.
smallThresh = cMax(J) / maxRatio;
removeEntry = ismember(J, badCols) & (absV < smallThresh);

%Keep all entries that are not removed
keep = ~removeEntry;
newI = I(keep);
newJ = J(keep);
newV = V(keep);

%Auxiliary entries are accumulated here
addI = zeros(0,1); addJ = zeros(0,1); addV = zeros(0,1);
convLb = zeros(0,1); convUb = zeros(0,1);
nStar = 0;   %auxiliary metabolite rows added so far
nConv = 0;   %conversion variable columns added so far

remIdx = find(removeEntry);
for k = 1:numel(remIdx)
    e         = remIdx(k);
    metRow    = I(e);
    parentCol = J(e);
    coef      = V(e);
    curCmax   = cMax(parentCol);
    %Deliver the (signed) coefficient "coef" of metabolite "metRow" into
    %column "parentCol", splitting recursively until it is within range
    while true
        small = curCmax / maxRatio;
        if abs(coef) >= small
            %Coefficient is acceptable; place it directly
            addI(end+1,1) = metRow;    addJ(end+1,1) = parentCol; addV(end+1,1) = coef; %#ok<AGROW>
            break
        end
        %Otherwise introduce a star metabolite and a conversion variable
        beta    = sign(coef) * small;
        nStar   = nStar + 1; starRow = nOrigRows + nStar;
        nConv   = nConv + 1; convCol = nOrigCols + nConv;
        %The parent column now consumes the star metabolite (coefficient beta)
        addI(end+1,1) = starRow; addJ(end+1,1) = parentCol; addV(end+1,1) = beta; %#ok<AGROW>
        %The conversion variable produces one unit of the star metabolite
        addI(end+1,1) = starRow; addJ(end+1,1) = convCol;   addV(end+1,1) = 1;    %#ok<AGROW>
        %Conversion variables are determined by the equality and left unbounded
        convLb(end+1,1) = -inf; convUb(end+1,1) = inf; %#ok<AGROW>
        %Deliver the residual coefficient into the new conversion column. Its
        %reference scale is 1 (the production coefficient of the star met).
        coef      = -coef / beta;
        parentCol = convCol;
        curCmax   = 1;
    end
end

condInfo.applied = true;
condInfo.nSplits = nStar;

%Assemble the augmented matrix and extend the problem vectors. New rows are
%appended at the bottom, new columns at the right.
prob.A = sparse([newI; addI], [newJ; addJ], [newV; addV], ...
    nOrigRows + nStar, nOrigCols + nConv);
prob.b       = [prob.b(:);       zeros(nStar,1)];
prob.csense  = [prob.csense(:).', repmat('E',1,nStar)];
prob.c       = [prob.c(:);       zeros(nConv,1)];
prob.lb      = [prob.lb(:);      convLb];
prob.ub      = [prob.ub(:);      convUb];
prob.vartype = [prob.vartype(:).', repmat('C',1,nConv)];
end
