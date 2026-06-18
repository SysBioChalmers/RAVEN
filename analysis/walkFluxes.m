function walkFluxes(model, fluxes, startRxn, varargin)
% walkFluxes  Interactively navigate a flux distribution reaction by reaction.
%
% Starting from a chosen reaction, displays all flux-carrying reactions
% connected through shared metabolites, grouped by metabolite and labelled
% with each neighbour's role (produces / consumes). Select a neighbour by
% number to step to it; build a history and rewind with 'b'.
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure.
% fluxes : double
%     flux vector (length = numel(model.rxns)).
% startRxn : char or double
%     starting reaction ID (string) or column index.
%
% Name-Value Arguments
% --------------------
% cutoff : double
%     minimum |flux| to display a reaction as a neighbour (default 1e-8).
% maxPerMet : double
%     maximum neighbour reactions displayed per metabolite (default 8).
%
% Examples
% --------
%     sol = solveLP(model);
%     walkFluxes(model, sol.x, 'BIOMASS_Ec_iJO1366_core_53p95M')
%     walkFluxes(model, sol.x, 'PFK', 'cutoff', 1e-6)
%
% Notes
% -----
% Navigation commands:
%   <number>  — step to that neighbour reaction
%   b         — go back to the previous reaction
%   q         — quit

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'maxPerMet',8});
cutoff    = p.cutoff;
maxPerMet = p.maxPerMet;

rxnIdx  = resolveRxn_(model, startRxn);
history = zeros(0,1);
W       = 68;

while true
    % ---- Header ----
    fprintf('\n%s\n', repmat('=', 1, W));
    f   = fluxes(rxnIdx);
    rID = model.rxns{rxnIdx};
    rName = getRxnName_(model, rxnIdx);
    if isempty(rName)
        fprintf('  [%s]  flux: %+.6g\n', rID, f);
    else
        fprintf('  [%s]  %s\n  flux: %+.6g\n', rID, rName, f);
    end
    eqn = constructEquations(model, model.rxns(rxnIdx));
    fprintf('  %s\n', eqn{1});
    if isfield(model,'grRules') && rxnIdx <= numel(model.grRules) ...
            && ~isempty(model.grRules{rxnIdx})
        fprintf('  genes: %s\n', model.grRules{rxnIdx});
    end
    fprintf('%s\n', repmat('=', 1, W));

    % ---- Collect neighbours grouped by metabolite ----
    col         = full(model.S(:, rxnIdx));
    nonzeroMets = find(abs(col) > 0);

    % One number per unique reaction across all metabolite groups
    seenRxns   = containers.Map('KeyType','double','ValueType','double');
    allNeighbors = [];

    for mi = 1:numel(nonzeroMets)
        m    = nonzeroMets(mi);
        netF = col(m) * f;   % >0: current rxn produces m; <0: consumes m

        if abs(netF) < cutoff, continue; end

        if netF < 0
            role = 'consumed';
        else
            role = 'produced';
        end

        % Other flux-carrying reactions involving this metabolite
        row   = full(model.S(m, :));
        carry = find(abs(row .* fluxes') > cutoff);
        carry = carry(carry ~= rxnIdx);
        if isempty(carry), continue; end

        [~, si] = sort(abs(fluxes(carry)), 'descend');
        carry   = carry(si(1:min(end, maxPerMet)));

        mLbl = getMetLabel_(model, m);
        fprintf('\n  %s  [%s, %.4g]\n', mLbl, role, abs(netF));

        for ci = 1:numel(carry)
            nr         = carry(ci);
            neighNetF  = row(nr) * fluxes(nr);  % >0: nr produces m; <0: consumes m
            neighRole  = 'produces';
            if neighNetF < 0, neighRole = 'consumes'; end

            nrName = getRxnName_(model, nr);
            if numel(nrName) > 28, nrName = [nrName(1:25) '...']; end

            if isKey(seenRxns, nr)
                n = seenRxns(nr);   % re-use existing number; don't re-add
            else
                n = numel(allNeighbors) + 1;
                seenRxns(nr)     = n;
                allNeighbors(end+1) = nr; %#ok<AGROW>
            end

            fprintf('  %3d. %-14s  %+9.4g  %-8s  %s\n', ...
                n, model.rxns{nr}, fluxes(nr), neighRole, nrName);
        end
    end

    fprintf('\n');
    nN = numel(allNeighbors);

    if nN == 0
        fprintf('  (no flux-carrying neighbours at this cutoff)\n');
        fprintf('  b: go back   q: quit\n');
        resp = strtrim(input('  > ', 's'));
    else
        fprintf('  1-%d: step to reaction   b: go back   q: quit\n', nN);
        resp = strtrim(input('  > ', 's'));
    end

    if strcmpi(resp, 'q')
        fprintf('Navigator closed.\n');
        return

    elseif strcmpi(resp, 'b')
        if isempty(history)
            fprintf('  (no history)\n');
        else
            rxnIdx     = history(end);
            history(end) = [];
        end

    else
        n = str2double(resp);
        if ~isnan(n) && n >= 1 && n <= nN && n == round(n)
            history(end+1) = rxnIdx; %#ok<AGROW>
            rxnIdx         = allNeighbors(n);
        else
            fprintf('  Enter a number 1-%d, ''b'', or ''q''.\n', max(nN,1));
        end
    end
end
end

%--------------------------------------------------------------------------
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

function s = getRxnName_(model, rxnIdx)
    s = '';
    if isfield(model,'rxnNames') && rxnIdx <= numel(model.rxnNames)
        s = char(model.rxnNames{rxnIdx});
    end
end

function s = getMetLabel_(model, metIdx)
    if isfield(model,'metNames') && ~isempty(model.metNames{metIdx})
        s = model.metNames{metIdx};
    else
        s = model.mets{metIdx};
    end
end
