%run this test case with the command
%results = runtests('conditioningTests.m')
function tests = conditioningTests
tests = functiontests(localfunctions);
end

function singleSplitTest(testCase)
%A column (reaction) with coefficients 10, 1 and 1e-8 has a ratio of 1e9,
%well above the threshold, and should be split exactly once.
maxRatio = 1e6;
prob = makeProb(sparse([10, 2, 0; 1, 0, 3; 1e-8, 1, 0]));

[probNew, condInfo] = splitProbForConditioning(prob, maxRatio);

verifyTrue(testCase, condInfo.applied)
verifyEqual(testCase, condInfo.nSplits, 1)
verifyEqual(testCase, condInfo.nOrigRows, 3)
verifyEqual(testCase, condInfo.nOrigCols, 3)

%No column may exceed maxRatio after the split
verifyLessThanOrEqual(testCase, maxColRatio(probNew.A), maxRatio*1.0001)

%The feasible region must be preserved and the problem vectors extended
verifyEquivalent(testCase, prob, probNew, condInfo)
end

function recursiveSplitTest(testCase)
%A coefficient of 1e-13 against a max of 10 has a ratio of 1e12, which needs
%two levels of splitting when maxRatio is 1e6.
maxRatio = 1e6;
prob = makeProb(sparse([10, 0; 1e-13, 4]));

[probNew, condInfo] = splitProbForConditioning(prob, maxRatio);

verifyTrue(testCase, condInfo.applied)
verifyGreaterThanOrEqual(testCase, condInfo.nSplits, 2)
verifyLessThanOrEqual(testCase, maxColRatio(probNew.A), maxRatio*1.0001)
verifyEquivalent(testCase, prob, probNew, condInfo)
end

function noSplitTest(testCase)
%A well-scaled problem must be returned untouched.
maxRatio = 1e6;
prob = makeProb(sparse([10, 2; 1, 3]));

[probNew, condInfo] = splitProbForConditioning(prob, maxRatio);

verifyFalse(testCase, condInfo.applied)
verifyEqual(testCase, probNew.A, prob.A)
end

function idempotenceTest(testCase)
%Splitting an already-conditioned problem a second time changes nothing.
maxRatio = 1e6;
prob = makeProb(sparse([10, 2, 0; 1, 0, 3; 1e-8, 1, 0]));

[probNew, ~]      = splitProbForConditioning(prob, maxRatio);
[~, condInfo2]    = splitProbForConditioning(probNew, maxRatio);

verifyFalse(testCase, condInfo2.applied)
end

% --- helpers --------------------------------------------------------------

function prob = makeProb(A)
[m,n]        = size(A);
prob.A       = A;
prob.b       = zeros(m,1);
prob.c       = [1; zeros(n-1,1)];
prob.lb      = -1000*ones(n,1);
prob.ub      =  1000*ones(n,1);
prob.csense  = repmat('E',1,m);
prob.osense  = 1;
prob.vartype = repmat('C',1,n);
end

function r = maxColRatio(A)
A = abs(A);
r = 0;
for j = 1:size(A,2)
    v = nonzeros(A(:,j));
    if numel(v) >= 2
        r = max(r, max(v)/min(v));
    end
end
end

function verifyEquivalent(testCase, prob, probNew, condInfo)
%Verify that the augmented problem has the same feasible region in the
%original variables, and that all problem vectors were extended consistently.
nR = condInfo.nOrigRows;
nC = condInfo.nOrigCols;
nS = condInfo.nSplits;

Aaug = full(probNew.A);
TL = Aaug(1:nR,      1:nC);
TR = Aaug(1:nR,      nC+1:end);
BL = Aaug(nR+1:end,  1:nC);
BR = Aaug(nR+1:end,  nC+1:end);

%The conversion variables are pinned by the star-metabolite rows:
%  BL*x + BR*conv = 0  ->  conv = -BR\BL * x
%so the original constraint rows evaluate to (TL + TR*(-BR\BL))*x, which must
%reproduce the original constraint matrix exactly.
M = -BR\BL;
verifyLessThan(testCase, max(abs((TL + TR*M) - full(prob.A)),[],'all'), 1e-9)

verifyEqual(testCase, numel(probNew.b),       nR+nS)
verifyEqual(testCase, numel(probNew.csense),  nR+nS)
verifyEqual(testCase, numel(probNew.c),       nC+nS)
verifyEqual(testCase, numel(probNew.lb),      nC+nS)
verifyEqual(testCase, numel(probNew.ub),      nC+nS)
verifyEqual(testCase, numel(probNew.vartype), nC+nS)
end
