function [pathRxns, pathMets, cumFrac] = traceFluxPath(model, fluxes, fromRxn, toRxn, varargin)
% traceFluxPath  Find the highest-flux-fraction path between two reactions.
%
% Traces how flux flows from one reaction to another through the metabolic
% network, quantifying what fraction of the source reaction's flux reaches
% the target. At each metabolite junction the flux is split proportionally
% among all consuming reactions; the path with the highest cumulative
% fraction is returned and printed.
%
% By default the search excludes currency metabolites (ATP, NAD, H2O, etc.)
% so the returned path represents material (carbon-skeleton) flow rather
% than cofactor shortcuts. Set traceMaterial=false to disable this.
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
% traceMaterial : logical
%     exclude common currency/cofactor metabolites (ATP, ADP, NAD, NADH,
%     H2O, CoA, Pi, CO2, …) so the path follows material (carbon-skeleton)
%     flux rather than energy-carrier shortcuts (default true).
% carbonOnly : logical
%     additionally require each intermediate metabolite to contain at least
%     one carbon atom, inferred from model.metFormulas. Drops H+, H2O, Pi,
%     O2, NH3, CO2 that are not already caught by the currency list.
%     Silently ignored when model.metFormulas is absent (default false).
% excludeMets : cell of char
%     extra metabolite names or IDs (after stripping compartment suffix) to
%     exclude, in addition to the built-in currency list (default {}).
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
%     % disable material filter to see all connections including ATP/ADP:
%     [p, m, f] = traceFluxPath(model, sol.x, 'PFK', 'ATPS4rpp', 'traceMaterial', false)
%     % add model-specific cofactors to the exclusion list:
%     [p, m, f] = traceFluxPath(model, sol.x, 'PFK', 'CS', 'excludeMets', {'acetyl-CoA'})
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

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'maxHops',30; 'verbose',true; ...
    'traceMaterial',true; 'carbonOnly',false; 'excludeMets',{}});
cutoff        = p.cutoff;
maxHops       = p.maxHops;
verbose       = p.verbose;
traceMaterial = p.traceMaterial;
carbonOnly    = p.carbonOnly;
excludeMets   = p.excludeMets;
if ischar(excludeMets), excludeMets = {excludeMets}; end

% ---- Build exclusion set ----
if traceMaterial
    currencyList = [defaultCurrencyMets_(), excludeMets(:)'];
else
    currencyList = excludeMets(:)';
end

% Pre-compute carbon counts per metabolite (for carbonOnly filter)
hasCarbonFilter = false;
carbonCounts    = [];
if carbonOnly
    if isfield(model,'metFormulas') && ~isempty(model.metFormulas)
        carbonCounts    = cellfun(@countCarbon_, model.metFormulas);
        hasCarbonFilter = true;
    else
        warning('traceFluxPath:noFormulas', ...
            'carbonOnly requires model.metFormulas; the option is ignored.');
    end
end

% ---- Resolve reactions ----
pathRxns = {};
pathMets = {};
cumFrac  = 0;

fromIdx = resolveRxn_(model, fromRxn);
toIdx   = resolveRxn_(model, toRxn);

if fromIdx == toIdx
    pathRxns = model.rxns(fromIdx);
    pathMets = {};
    cumFrac  = 1;
    if verbose, printPath_(model, fluxes, pathRxns, pathMets, cumFrac, ...
            fromRxn, toRxn, traceMaterial || ~isempty(excludeMets) || carbonOnly); end
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
    col  = full(model.S(:, cur));
    fcur = fluxes(cur);
    for m = 1:numel(col)
        if abs(col(m)) < 1e-15, continue; end
        if col(m) * fcur <= 0, continue; end   % skip consumed or zero-contribution mets

        % ---- Material-flux filters ----
        if ~isempty(currencyList) && isCurrencyMet_(model, m, currencyList)
            continue
        end
        if hasCarbonFilter && carbonCounts(m) == 0
            continue
        end

        % Reactions that net-consume this metabolite
        row  = full(model.S(m, :));
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
    materialMode = traceMaterial || ~isempty(excludeMets) || carbonOnly;
    printPath_(model, fluxes, pathRxns, pathMets, cumFrac, fromRxn, toRxn, materialMode);
end
end

%--------------------------------------------------------------------------
function e = makeEntry_(rxn, frac, rpath, mpath)
    e.rxn   = rxn;
    e.frac  = frac;
    e.rpath = rpath;
    e.mpath = mpath;
end

function printPath_(model, fluxes, pathRxns, pathMets, cumFrac, fromRxn, toRxn, materialMode)
    W   = 66;
    bar = repmat('=', 1, W);
    fprintf('\n%s\n', bar);
    fprintf('  Flux path: %s \x2192 %s\n', fromRxn, toRxn);
    if materialMode
        fprintf('  (material flux — currency metabolites excluded)\n');
    end
    fprintf('%s\n', bar);

    if isempty(pathRxns)
        fprintf('  No forward flux path found between these reactions.\n');
        if materialMode
            fprintf('  Tip: try traceMaterial=false to include cofactor connections,\n');
            fprintf('       or excludeMets={} to see what the filter removed.\n');
        else
            fprintf('  (They may be in parallel branches or connected only in reverse.)\n');
        end
        fprintf('%s\n\n', bar);
        return
    end

    % One-line path with junction fractions
    parts = pathRxns{1};
    for i = 2:numel(pathRxns)
        mi = find(strcmp(model.mets, pathMets{i-1}));
        row = full(model.S(mi, :));
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
        ri  = find(strcmp(model.rxns, pathRxns{i}));
        nm  = '';
        if isfield(model,'rxnNames') && ri <= numel(model.rxnNames)
            nm = char(model.rxnNames{ri});
        end
        if numel(nm) > 25, nm = [nm(1:22) '...']; end
        fprintf('  %2d. %-14s  flux: %+10.4g   %-25s  %s\n', i, pathRxns{i}, fluxes(ri), nm, eqns{i});
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

function tf = isCurrencyMet_(model, m, currencyList)
% Returns true if metabolite m matches any entry in currencyList.
% Matches against metNames (case-insensitive exact) and against the
% metabolite ID after stripping compartment suffixes like [c], [m], _c.
    targets = {};
    if isfield(model,'metNames') && m <= numel(model.metNames) && ~isempty(model.metNames{m})
        targets{end+1} = lower(strtrim(model.metNames{m}));
    end
    % Strip compartment from ID: [c], (c), _c, _C at end
    rawId = model.mets{m};
    stripped = regexprep(rawId, '[\[\(][^\]\)]*[\]\)]$|_[a-zA-Z]\d*$', '');
    targets{end+1} = lower(stripped);

    for i = 1:numel(currencyList)
        p = lower(strtrim(currencyList{i}));
        for j = 1:numel(targets)
            if strcmp(targets{j}, p)
                tf = true;
                return
            end
        end
    end
    tf = false;
end

function nC = countCarbon_(formula)
% Count carbon atoms in a molecular formula string like 'C6H12O6'.
    if isempty(formula)
        nC = 0; return
    end
    tok = regexp(char(formula), 'C(\d*)', 'tokens', 'once');
    if isempty(tok)
        nC = 0;
    elseif isempty(tok{1})
        nC = 1;
    else
        nC = str2double(tok{1});
    end
end

function lst = defaultCurrencyMets_()
% Common currency / cofactor metabolites matched case-insensitively.
% Bare comment lines inside {} act as row separators in older MATLAB, so
% all comments are placed after ... on data lines instead.
    lst = { ...
        'ATP','ADP','AMP','dATP','dADP','dAMP', ... % adenosine phosphates
        'atp','adp','amp','datp','dadp','damp', ...
        'GTP','GDP','GMP','dGTP','dGDP','dGMP', ... % other NTPs/NDPs/NMPs
        'CTP','CDP','CMP','dCTP','dCDP','dCMP', ...
        'UTP','UDP','UMP','dUTP','dUDP','dUMP', ...
        'TTP','TDP','TMP','dTTP','dTDP','dTMP', ...
        'gtp','gdp','gmp','dgtp','dgdp','dgmp', ...
        'ctp','cdp','cmp','dctp','dcdp','dcmp', ...
        'utp','udp','ump','dutp','dudp','dump', ...
        'ttp','tdp','tmp','dttp','dtdp','dtmp', ...
        'NAD','NADH','NADP','NADPH','NAD+','NADP+', ... % nicotinamide cofactors
        'nad','nadh','nadp','nadph', ...
        'FAD','FADH2','FMN','FMNH2', ... % flavin cofactors
        'fad','fadh2','fmn','fmnh2', ...
        'CoA','coenzyme A','coenzyme a','Coenzyme A', ... % coenzyme A
        'coa', ...
        'H2O','water','H+','proton','OH-','hydroxide', ... % water / proton
        'h2o','h','oh1', ...
        'phosphate','orthophosphate','pyrophosphate','diphosphate', ... % inorganic phosphate
        'pi','ppi','pp', ...
        'CO2','carbon dioxide','bicarbonate','HCO3-', ... % carbon dioxide
        'co2','hco3', ...
        'O2','oxygen','NH3','ammonia','NH4+','ammonium','nitrogen', ... % oxygen / nitrogen
        'o2','nh3','nh4', ...
        'sulfate','sulfite','sulfide','thiosulfate', ... % sulfur species
        'thioredoxin','thioredoxin-SH','thioredoxin-S2', ... % redox carriers
        'trdrd','trdox', ...
        'ferredoxin','ferredoxin reduced','ferredoxin oxidized', ...
        'fdred','fdox', ...
        'ubiquinol','ubiquinone','menaquinol','menaquinone', ...
        'q8h2','q8','mqn8','mql8', ...
        'lipoamide','dihydrolipoamide', ...
        'S-adenosyl-L-methionine','S-adenosylhomocysteine', ... % one-carbon donors
        'sam','sah','amet','ahcys', ...
        '5-methyltetrahydrofolate','tetrahydrofolate','dihydrofolate', ...
        'thf','dhf','mlthf','methf' ...
    };
end
