function [pathRxns, pathMets, cumFrac] = traceFluxPath(model, fluxes, fromRxn, toRxn, varargin)
% traceFluxPath  Find the highest-flux-fraction path between two reactions.
%
% Traces how flux flows from one reaction to another through the metabolic
% network, quantifying what fraction of the source reaction's flux reaches
% the target. At each metabolite junction the flux is split proportionally
% among all consuming reactions; the path with the highest cumulative
% fraction is returned and printed.
%
% This answers the question: "Which intermediate reactions are most
% responsible for routing flux from reaction A to reaction B?"
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure.
% fluxes : double
%     flux vector (length = numel(model.rxns)).
% fromRxn : char
%     source reaction ID.
% toRxn : char
%     target reaction ID.
%
% Name-Value Arguments
% --------------------
% cutoff : double
%     minimum |flux| for a reaction to be considered flux-carrying
%     (default 1e-8).
% maxHops : double
%     maximum number of reactions in the path; longer paths are pruned
%     (default 30).
% verbose : logical
%     print the path to the command window (default true).
%
% Returns
% -------
% pathRxns : cell
%     cell array of reaction IDs on the best path, including fromRxn and
%     toRxn. Empty if no path is found.
% pathMets : cell
%     metabolite IDs at each junction (numel = numel(pathRxns)-1).
% cumFrac : double
%     fraction of fromRxn's flux that routes through this exact path
%     (0 = no path found, 1 = all flux reaches toRxn directly).
%
% Examples
% --------
%     [p, m, f] = traceFluxPath(model, sol.x, 'PFK', 'ATPS4rpp')
%
% Notes
% -----
% The algorithm follows metabolites in the forward direction only: it
% traces the metabolites that fromRxn PRODUCES and finds which consuming
% reactions eventually reach toRxn. If fromRxn and toRxn are not connected
% in the forward flux direction (e.g. they are in parallel branches), an
% empty path is returned. A best-first (greedy) search ensures the
% highest-fraction path is found efficiently; a visited set prevents
% revisiting reactions.

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'maxHops',30; 'verbose',true});
cutoff  = p.cutoff;
maxHops = p.maxHops;
verbose = p.verbose;

pathRxns = {};
pathMets = {};
cumFrac  = 0;

fromIdx = resolveRxn_(model, fromRxn);
toIdx   = resolveRxn_(model, toRxn);

if fromIdx == toIdx
    pathRxns = model.rxns(fromIdx);
    pathMets = {};
    cumFrac  = 1;
    if verbose, printPath_(model, fluxes, pathRxns, pathMets, cumFrac, fromRxn, toRxn); end
    return
end

nRxns = numel(model.rxns);

% ---- Best-first search (highest cumulative fraction first) ----
% Queue entries are structs: .rxn (idx), .frac (cumFrac), .rpath, .mpath
bestAtNode = zeros(nRxns,1);
bestAtNode(fromIdx) = 1.0;

queue = makeEntry_(fromIdx, 1.0, fromIdx, zeros(0,1));

bestFrac  = 0;
bestRPath = [];
bestMPath = [];

while ~isempty(queue)
    % Pop highest-fraction entry
    entry = queue(1);
    queue(1) = [];

    cur   = entry.rxn;
    frac  = entry.frac;
    rpath = entry.rpath;
    mpath = entry.mpath;

    if cur == toIdx
        if frac > bestFrac
            bestFrac  = frac;
            bestRPath = rpath;
            bestMPath = mpath;
        end
        continue
    end

    if numel(rpath) > maxHops, continue; end

    % Find metabolites net-produced by cur
    col  = model.S(:, cur);
    fcur = fluxes(cur);
    for m = 1:numel(col)
        if abs(col(m)) < 1e-15, continue; end
        if col(m) * fcur <= 0, continue; end   % skip consumed or zero-contribution mets

        % Reactions that net-consume this metabolite
        row  = model.S(m, :);
        fVec = fluxes(:)';
        consumers = find(row .* fVec < -cutoff);
        consumers = consumers(~ismember(consumers, rpath));  % no revisiting
        if isempty(consumers), continue; end

        totalCons = sum(abs(row(consumers) .* fVec(consumers)));

        for ci = 1:numel(consumers)
            nr      = consumers(ci);
            frac_r  = abs(row(nr) * fluxes(nr)) / totalCons;
            newFrac = frac * frac_r;

            if newFrac > bestAtNode(nr) || nr == toIdx
                bestAtNode(nr) = max(bestAtNode(nr), newFrac);
                newEntry = makeEntry_(nr, newFrac, [rpath, nr], [mpath; m]);

                % Insert in decreasing-fraction order
                if isempty(queue)
                    queue = newEntry;
                else
                    fracs = [queue.frac];
                    pos   = find(fracs < newFrac, 1);
                    if isempty(pos)
                        queue(end+1) = newEntry; %#ok<AGROW>
                    else
                        queue = [queue(1:pos-1), newEntry, queue(pos:end)];
                    end
                end
            end
        end
    end
end

% ---- Package and display ----
cumFrac = bestFrac;
if ~isempty(bestRPath)
    pathRxns = model.rxns(bestRPath);
    pathMets = model.mets(bestMPath);
end

if verbose
    printPath_(model, fluxes, pathRxns, pathMets, cumFrac, fromRxn, toRxn);
end
end

%--------------------------------------------------------------------------
function e = makeEntry_(rxn, frac, rpath, mpath)
    e.rxn   = rxn;
    e.frac  = frac;
    e.rpath = rpath;
    e.mpath = mpath;
end

function printPath_(model, fluxes, pathRxns, pathMets, cumFrac, fromRxn, toRxn)
    W   = 66;
    bar = repmat('=', 1, W);
    fprintf('\n%s\n', bar);
    fprintf('  Flux path: %s \x2192 %s\n', fromRxn, toRxn);
    fprintf('%s\n', bar);

    if isempty(pathRxns)
        fprintf('  No forward flux path found between these reactions.\n');
        fprintf('  (They may be in parallel branches or connected only in reverse.)\n');
        fprintf('%s\n\n', bar);
        return
    end

    % One-line path with junction fractions
    parts = pathRxns{1};
    for i = 2:numel(pathRxns)
        mi = find(strcmp(model.mets, pathMets{i-1}));
        row = model.S(mi, :);
        consumers = find(row .* fluxes' < -1e-12);
        if isempty(consumers)
            pct = NaN;
        else
            totalCons = sum(abs(row(consumers) .* fluxes(consumers)'));
            ni = find(strcmp(model.rxns, pathRxns{i}));
            pct = abs(row(ni) * fluxes(ni)) / totalCons * 100;
        end
        mLbl = getMetLabel_(model, mi);
        if isnan(pct)
            parts = [parts sprintf(' \x2192[%s]\x2192 %s', mLbl, pathRxns{i})]; %#ok<AGROW>
        else
            parts = [parts sprintf(' \x2192[%s: %.0f%%]\x2192 %s', mLbl, pct, pathRxns{i})]; %#ok<AGROW>
        end
    end
    fprintf('\n  %s\n\n', parts);

    fprintf('  Cumulative fraction: %.2f%%\n', cumFrac*100);
    fprintf('  Path length:         %d reactions\n', numel(pathRxns));

    fprintf('\n  Reactions:\n');
    eqns = constructEquations(model, pathRxns);
    for i = 1:numel(pathRxns)
        ri = find(strcmp(model.rxns, pathRxns{i}));
        fprintf('  %2d. %-14s  flux: %+10.4g   %s\n', i, pathRxns{i}, fluxes(ri), eqns{i});
    end
    fprintf('\n%s\n\n', bar);
end

function idx = resolveRxn_(model, rxn)
    if isnumeric(rxn)
        idx = rxn;
    else
        idx = find(strcmp(model.rxns, char(rxn)));
        if isempty(idx)
            error('Reaction ''%s'' not found in model.', rxn);
        end
        idx = idx(1);
    end
end

function s = getMetLabel_(model, metIdx)
    if isfield(model,'metNames') && ~isempty(model.metNames{metIdx})
        s = model.metNames{metIdx};
    else
        s = model.mets{metIdx};
    end
end
