function [x, E, converged] = sampleMaxVolEllipse(A, b, x0, maxiter, tol, reg)
% sampleMaxVolEllipse  Maximum-volume inscribed ellipsoid of a polytope.
%
% Solves  max log det E  over centre x and SPD matrix E such that the
% ellipsoid {x + E*s : ||s||_2 <= 1} is contained in the polytope
% {z : A*z <= b}, using the Zhang & Gao (2003) primal-dual interior-point
% method (the regularised variant used in COBRA's CHRR sampler). This is the
% rounding step of sampleCHRR: the returned transform E maps the polytope to a
% well-conditioned body in which coordinate hit-and-run mixes rapidly.
%
% This function is pure linear algebra (no LP solver) except when x0 is empty,
% in which case a Chebyshev centre is obtained via sampleChebyshevCenter.
%
% Parameters
% ----------
% A : double
%     m-by-n constraint matrix of the polytope {z : A*z <= b}. Requires
%     m >= n+1 and a non-empty interior.
% b : double
%     m-by-1 right-hand side.
% x0 : double (optional)
%     a strictly interior point (A*x0 < b). If empty, computed by
%     sampleChebyshevCenter(A, b).
% maxiter : double (optional)
%     maximum interior-point iterations (default 150).
% tol : double (optional)
%     convergence tolerance (default 1e-6).
% reg : double (optional)
%     diagonal / Levenberg regularisation keeping the Newton solves
%     well-conditioned (default 1e-8).
%
% Returns
% -------
% x : double
%     n-by-1 ellipsoid centre.
% E : double
%     n-by-n lower-triangular rounding transform with E*E' = E2 (the SPD
%     ellipsoid matrix). The rounding substitution is z = x + E*y.
% converged : logical
%     whether tol was reached within maxiter. A false value is not fatal —
%     the last iterate is still a valid (if looser) rounding.
%
% Notes
% -----
% Validated against analytic cases: a box maps to the unit ball (E = I); a
% sheared box {z : -1 <= M*z <= 1} gives E*E' = inv(M)*inv(M)'; a triangle
% gives its Steiner inellipse.
%
% Based on Y. Zhang & L. Gao (2003), SIAM J. Optim. 14:53-76; as used in
% Haraldsdottir et al. (2017) Bioinformatics 33:1741 (CHRR).
%
% See also
% --------
% sampleCHRR, sampleChebyshevCenter, randomSampling

if nargin < 4 || isempty(maxiter); maxiter = 150; end
if nargin < 5 || isempty(tol);     tol = 1e-6;     end
if nargin < 6 || isempty(reg);     reg = 1e-8;     end

[m, n] = size(A);
b = b(:);
% bnrm is the norm of the ORIGINAL b (before the row-normalisation below sets
% b <- ones). This matches COBRA's chrrSampler reference; do not "correct" it to
% sqrt(m), or convergence scaling diverges from the reference (and from the
% validated Python max_volume_ellipsoid this file mirrors).
bnrm = norm(b);
minmu = 1e-8;
tau0  = 0.75;

if nargin < 3 || isempty(x0)
    x0 = sampleChebyshevCenter(A, b);
end
x0 = x0(:);

bmAx0 = b - A*x0;
if any(bmAx0 <= 0)
    error('RAVEN:sampling', 'sampleMaxVolEllipse: x0 is not strictly interior (A*x0 < b).');
end

% Row-normalise so the constraint RHS is all ones; solve in shifted coords.
A = A ./ bmAx0;
b = ones(m, 1);

x = zeros(n, 1);
y = ones(m, 1);
bmAx = b;
z = zeros(m, 1);
astep = 0;
Adx = zeros(n, 1);
converged = false;
E2 = eye(n);

for it = 1:maxiter
    if it > 1
        bmAx = bmAx - astep*Adx;
    end

    AtYA = A' * (A .* y);          % A' * diag(y) * A   (n-by-n)
    E2 = inv(AtYA);                %#ok<MINV> explicit inverse mirrors the reference
    Q  = A * E2 * A';              % m-by-m
    h  = sqrt(max(diag(Q), 0));

    if it == 1
        t = min(bmAx ./ h);
        y = y / t^2;
        h = t * h;
        z = max(1e-1, bmAx - h);
        Q = t^2 * Q;
    end

    yz = y .* z;
    yh = y .* h;
    gap = sum(yz) / m;
    rmu = max(min(0.5, gap)*gap, minmu);

    R1 = -A' * yh;                 % dual residual      (n-by-1)
    R2 = bmAx - h - z;             % primal/slack       (m-by-1)
    R3 = rmu - yz;                 % complementarity    (m-by-1)
    res = max([max(abs(R1)), max(abs(R2)), max(abs(R3))]);

    if res < tol*(1 + bnrm) && rmu <= minmu
        x = x + x0;
        converged = true;
        break;
    end

    YQ   = y .* Q;                 % diag(y) * Q
    YQQY = YQ .* YQ';              % Hadamard product
    y2h  = 2 * yh;
    G    = YQQY + diag(max(reg, y2h .* z));
    YA   = y .* A;                 % diag(y) * A
    T    = G \ ((h + z) .* YA);    % m-by-n
    ATP  = (y2h .* T - YA)';       % n-by-m
    R3Dy = R3 ./ y;
    R23  = R2 - R3Dy;
    ATP_A = ATP*A + reg*eye(n);
    dx   = ATP_A \ (R1 + ATP*R23);
    Adx  = A * dx;
    dyDy = G \ (y2h .* (Adx - R23));
    dy   = y .* dyDy;
    dz   = R3Dy - z .* dyDy;

    ax = -1 / min([-Adx ./ bmAx; -0.5]);
    ay = -1 / min([ dyDy;        -0.5]);
    az = -1 / min([ dz ./ z;     -0.5]);
    tau = max(tau0, 1 - res);
    astep = tau * min([1, ax, ay, az]);

    x = x + astep*dx;
    y = y + astep*dy;
    z = z + astep*dz;
end

if ~converged
    x = x + x0;
end

E = nearestSPDchol(E2);
end


function E = nearestSPDchol(M)
% Lower-triangular Cholesky factor of the nearest SPD matrix to M, so E*E' = M.
M = (M + M') / 2;
[E, p] = chol(M, 'lower');
if p == 0
    return;
end
ev = eig(M);
jitter = max(-min(ev), 0) + 1e-12;
nn = size(M, 1);
for k = 1:30
    [E, p] = chol(M + jitter*eye(nn), 'lower');
    if p == 0
        return;
    end
    jitter = jitter * 10;
end
E = chol(M + jitter*eye(nn), 'lower');
end
