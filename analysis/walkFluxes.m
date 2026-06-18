function walkFluxes(model, fluxes, startRxn, varargin)
% walkFluxes  Interactively navigate a flux distribution reaction by reaction.
%
% Starting from a chosen reaction, displays all flux-carrying reactions
% connected through shared metabolites. The user selects a neighbour by
% number to step to that reaction, building a navigation history that can
% be rewound with 'b'. Use this to trace how flux flows through a network
% without knowing the full topology in advance.
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
% Navigation commands at the prompt:
%   <number>  — step to that neighbour reaction
%   b         — go back to the previous reaction
%   q         — quit

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'maxPerMet',8});
cutoff    = p.cutoff;
maxPerMet = p.maxPerMet;

% Resolve starting reaction
rxnIdx = resolveRxn_(model, startRxn);

history = zeros(0,1);

fprintf('\nFlux navigator — b=back, q=quit\n');

while true
    % ---- Header ----
    W  = 66;
    fprintf('\n%s\n', repmat('=', 1, W));
    f = fluxes(rxnIdx);
    rID = model.rxns{rxnIdx};
    rName = getRxnName_(model, rxnIdx);
    if isempty(rName)
        fprintf('  [%s]   flux: %+.6g\n', rID, f);
    else
        fprintf('  [%s]  %s\n  flux: %+.6g\n', rID, rName, f);
    end
    eqn = constructEquations(model, model.rxns(rxnIdx));
    fprintf('  %s\n', eqn{1});
    fprintf('%s\n', repmat('=', 1, W));

    % ---- Collect neighbours grouped by metabolite ----
    col = model.S(:, rxnIdx);
    nonzeroMets = find(col);

    % Track which reactions we've already assigned a number
    seenRxns = containers.Map('KeyType','double','ValueType','double');
    allNeighbors = [];   % reaction indices in display order (one entry per number)

    for mi = 1:numel(nonzeroMets)
        m    = nonzeroMets(mi);
        s_m  = col(m);
        netF = s_m * f;         % positive = produced here, negative = consumed here

        if netF < 0
            role = 'consumed';
        elseif netF > 0
            role = 'produced';
        else
            role = 'balanced';  % stoich but zero flux contribution
        end

        % Reactions involving this metabolite with non-trivial flux
        row   = model.S(m, :);
        carry = find(abs(row .* fluxes') > cutoff);
        carry = carry(carry ~= rxnIdx);
        if isempty(carry), continue; end

        [~, si] = sort(abs(fluxes(carry)), 'descend');
        carry = carry(si(1:min(end, maxPerMet)));

        % Print metabolite header
        mLbl = getMetLabel_(model, m);
        fprintf('\n  %s [%s here, net %+.4g]:\n', mLbl, role, netF);

        % Print each neighbour
        for ci = 1:numel(carry)
            nr = carry(ci);
            if isKey(seenRxns, nr)
                % Already numbered — show cross-reference
                existN = seenRxns(nr);
                fprintf('     %-12s  flux: %+10.4g   (same as #%d)\n', ...
                    model.rxns{nr}, fluxes(nr), existN);
            else
                n = numel(allNeighbors) + 1;
                seenRxns(nr) = n;
                allNeighbors(end+1) = nr; %#ok<AGROW>
                nrName = getRxnName_(model, nr);
                if numel(nrName) > 32, nrName = [nrName(1:29) '...']; end
                fprintf('  %3d. %-12s  flux: %+10.4g   %s\n', ...
                    n, model.rxns{nr}, fluxes(nr), nrName);
            end
        end
    end

    fprintf('\n');

    if isempty(allNeighbors)
        fprintf('  (no flux-carrying neighbours at this cutoff)\n');
        resp = strtrim(input('  b=back, q=quit: ', 's'));
    else
        resp = strtrim(input('  Navigate (number / b / q): ', 's'));
    end

    if strcmpi(resp, 'q')
        fprintf('Navigator closed.\n');
        return

    elseif strcmpi(resp, 'b')
        if isempty(history)
            fprintf('  No history.\n');
        else
            rxnIdx = history(end);
            history(end) = [];
        end

    else
        n = str2double(resp);
        if ~isnan(n) && n >= 1 && n <= numel(allNeighbors) && n == round(n)
            history(end+1) = rxnIdx; %#ok<AGROW>
            rxnIdx = allNeighbors(n);
        else
            if isempty(allNeighbors)
                fprintf('  Enter ''b'' or ''q''.\n');
            else
                fprintf('  Enter a number 1\x2013%d, ''b'', or ''q''.\n', numel(allNeighbors));
            end
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
