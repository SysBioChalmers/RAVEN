function result = compareFluxes(model, fluxes1, fluxes2, varargin)
% compareFluxes  Find the most significant flux changes between two conditions.
%
% Compares two flux vectors from the same model and returns all reactions
% that changed, sorted by absolute flux difference. Classifies reactions as
% turned on, turned off, reversed direction, or changed magnitude.
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure.
% fluxes1 : double
%     reference flux vector (condition A, e.g. wild-type).
% fluxes2 : double
%     comparison flux vector (condition B, e.g. mutant or knockin).
%
% Name-Value Arguments
% --------------------
% cutoff : double
%     minimum |flux| to consider a reaction active (default 1e-8).
% nMax : double
%     maximum rows printed in the summary table (default 20).
% verbose : logical
%     print the summary table (default true).
%
% Returns
% -------
% result : struct
%     .turnedOn  — cell, rxn IDs newly active in fluxes2
%     .turnedOff — cell, rxn IDs silenced in fluxes2
%     .flipped   — cell, rxn IDs that reversed sign between conditions
%     .changed   — struct with fields (all reactions with |dflux|>cutoff,
%                  sorted by |dflux| descending):
%                    .rxn       — cell of reaction IDs
%                    .flux1     — double, flux in condition A
%                    .flux2     — double, flux in condition B
%                    .absDelta  — double, |flux2 - flux1|
%                    .relChange — double, (flux2-flux1)/|flux1|; NaN when
%                                 flux1 is near-zero (reaction emerged)
%                    .type      — cell, 'on'|'off'|'flip'|'' per reaction
%
% Examples
% --------
%     solA = solveLP(model);
%     modelAna = setParam(model, 'ub', 'O2exchange', 0);
%     solB = solveLP(modelAna);
%     result = compareFluxes(model, solA.x, solB.x)
%     result = compareFluxes(model, solA.x, solB.x, 'nMax', 50)
%
% See also
%     walkFluxes, traceFluxPath, modelSummary

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'nMax',20; 'verbose',true});
cutoff  = p.cutoff;
nMax    = p.nMax;
verbose = p.verbose;

f1    = fluxes1(:);
f2    = fluxes2(:);
nRxns = numel(model.rxns);

if numel(f1) ~= nRxns || numel(f2) ~= nRxns
    error('RAVEN:badInput', ...
        'Flux vectors must have %d elements (numel(model.rxns)).', nRxns);
end

active1 = abs(f1) > cutoff;
active2 = abs(f2) > cutoff;
delta   = f2 - f1;

onMask   = ~active1 &  active2;
offMask  =  active1 & ~active2;
flipMask =  active1 &  active2 & sign(f1) ~= sign(f2);

result.turnedOn  = model.rxns(onMask);
result.turnedOff = model.rxns(offMask);
result.flipped   = model.rxns(flipMask);

% All changed reactions sorted by |delta| descending
changedMask = abs(delta) > cutoff;
[~, si]     = sort(abs(delta(changedMask)), 'descend');
idx = find(changedMask);
idx = idx(si);

% Relative change: (f2-f1)/|f1|; NaN when reference flux is near-zero
denom = abs(f1(idx));
denom(denom < cutoff) = NaN;
relChange = delta(idx) ./ denom;

% Type label per changed reaction
changedTypes = repmat({''}, numel(idx), 1);
for k = 1:numel(idx)
    ri = idx(k);
    if onMask(ri)
        changedTypes{k} = 'on';
    elseif offMask(ri)
        changedTypes{k} = 'off';
    elseif flipMask(ri)
        changedTypes{k} = 'flip';
    end
end

result.changed.rxn       = model.rxns(idx);
result.changed.flux1     = f1(idx);
result.changed.flux2     = f2(idx);
result.changed.absDelta  = abs(delta(idx));
result.changed.relChange = relChange;
result.changed.type      = changedTypes;

if verbose
    printTable_(model, result, nMax, nRxns);
end
end

%--------------------------------------------------------------------------
function printTable_(model, r, nMax, nTotal)
    W   = 72;
    bar = repmat('=', 1, W);
    sep = repmat('-', 1, W);
    nC  = numel(r.changed.rxn);

    fprintf('\n%s\n', bar);
    fprintf('  Flux comparison: %d of %d reactions changed\n', nC, nTotal);
    nOn  = numel(r.turnedOn);
    nOff = numel(r.turnedOff);
    nFlp = numel(r.flipped);
    if nOn + nOff + nFlp > 0
        fprintf('  (on: %d  |  off: %d  |  reversed: %d)\n', nOn, nOff, nFlp);
    end

    if nC == 0
        fprintf('%s\n\n', bar);
        return
    end

    nShow    = min(nC, nMax);
    showRxns = r.changed.rxn(1:nShow);
    eqns     = constructEquations(model, showRxns);

    fprintf('%s\n', sep);
    fprintf('  %-20s  %10s  %10s  %10s  %7s  %-5s  %-25s  %s\n', ...
        'Reaction', 'Cond.1', 'Cond.2', '|dflux|', 'drel%', 'Type', 'Name', 'Equation');
    fprintf('  %s\n', sep);

    for i = 1:nShow
        rel = r.changed.relChange(i);
        if isnan(rel)
            relStr = '    new';
        else
            relStr = sprintf('%+6.0f%%', rel * 100);
        end
        ri = find(strcmp(model.rxns, showRxns{i}), 1);
        nm = '';
        if isfield(model,'rxnNames') && ri <= numel(model.rxnNames)
            nm = char(model.rxnNames{ri});
        end
        if numel(nm) > 25, nm = [nm(1:22) '...']; end
        fprintf('  %-20s  %+10.4g  %+10.4g  %10.4g  %7s  %-5s  %-25s  %s\n', ...
            showRxns{i}, r.changed.flux1(i), r.changed.flux2(i), ...
            r.changed.absDelta(i), relStr, r.changed.type{i}, nm, eqns{i});
    end
    if nC > nMax
        fprintf('  ... (%d more, use nMax to show all)\n', nC - nMax);
    end
    fprintf('%s\n\n', bar);
end
