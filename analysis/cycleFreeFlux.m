function fluxOut = cycleFreeFlux(model, fluxes)
% cycleFreeFlux  Remove thermodynamically infeasible internal loops from a flux
%                distribution.
%
% Given a steady-state flux distribution, return an equivalent one that has the
% same exchange (boundary) fluxes but the smallest possible total internal flux
% and no net flux around internal loops. Method of Desouki et al. (2015),
% Bioinformatics 31:2159 (the LP-based loop removal cobra uses for its
% loopless='cycleFreeFlux' option): fix the boundary reactions and the sign of
% every internal reaction, then minimise the sum of absolute internal fluxes.
% A flux carried only to close a loop is removed because it adds to that sum
% without changing any exchange.
%
% Parameters
% ----------
% model : struct
%     a RAVEN model structure.
% fluxes : double
%     a flux vector (one value per reaction) to make cycle-free.
%
% Returns
% -------
% fluxOut : double
%     the loop-free flux vector, with the same boundary fluxes as the input.
%
% See also
% --------
% looplessFVA, getExchangeRxns, solveLP

tol = 1e-9;
nRxns = numel(model.rxns);
fluxes = fluxes(:);

m = model;
% Fix the boundary (exchange) reactions to their input values so the network's
% exchange with the environment is preserved.
[~, exchIdx] = getExchangeRxns(m);
isExch = false(nRxns,1); isExch(exchIdx) = true;
m.lb(exchIdx) = fluxes(exchIdx);
m.ub(exchIdx) = fluxes(exchIdx);

% Fix the sign of every internal reaction (so the cleaned solution keeps the
% same directions), and cap the magnitude at the input flux.
internal = find(~isExch);
for i = internal'
    v = fluxes(i);
    if v > tol
        m.lb(i) = max(m.lb(i), 0); m.ub(i) = min(m.ub(i), v);
    elseif v < -tol
        m.ub(i) = min(m.ub(i), 0); m.lb(i) = max(m.lb(i), v);
    else
        m.lb(i) = 0; m.ub(i) = 0;
    end
end

% Minimise the total internal flux (no objective; solveLP's minFlux pass
% minimises the sum of absolute fluxes, and the boundary fluxes are fixed).
m.c = zeros(nRxns,1);
sol = solveLP(m, 1);
if isempty(sol.x)
    fluxOut = fluxes;                 % infeasible under the fixings: return input
else
    fluxOut = sol.x;
end
end
