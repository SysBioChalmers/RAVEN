function x0 = sampleChebyshevCenter(A, b)
% sampleChebyshevCenter  Centre of the largest ball inscribed in a polytope.
%
% Solves the LP  max r  s.t.  a_i'*x + ||a_i||*r <= b_i,  r >= 0, giving a
% strictly interior point of {x : A*x <= b}. Used to initialise the
% interior-point MVE solver (sampleMaxVolEllipse) inside sampleCHRR.
%
% Parameters
% ----------
% A : double
%     m-by-n constraint matrix.
% b : double
%     m-by-1 right-hand side.
%
% Returns
% -------
% x0 : double
%     n-by-1 Chebyshev centre (strictly interior point).
%
% See also
% --------
% sampleMaxVolEllipse, sampleCHRR

[m, n] = size(A);
b = b(:);
norms = sqrt(sum(A.^2, 2));

% Variables [x (n, free); r (1, >=0)]; minimise -r (maximises the ball radius).
prob.A      = [A, norms];
prob.a      = prob.A;
prob.b      = b;
prob.csense = repmat('L', 1, m);
prob.c      = [zeros(n, 1); -1];
prob.osense = 1;                       % minimise
prob.lb     = [-inf(n, 1); 0];
prob.ub     = [ inf(n, 1); inf];
prob.vartype = repmat('C', 1, n + 1);

sol = optimizeProb(prob, [], false);
if ~checkSolution(sol)
    error('RAVEN:sampling', ...
        'sampleChebyshevCenter: LP infeasible — flux polytope has empty interior.');
end

x0 = sol.full(1:n);
r  = sol.full(n + 1);
if r <= 1e-12
    error('RAVEN:sampling', ...
        ['sampleChebyshevCenter: flux polytope has empty interior ' ...
         '(it is lower-dimensional than detected).']);
end
end
