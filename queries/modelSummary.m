function modelSummary(model, varargin)
% modelSummary  Print a concise overview of a RAVEN model and its flux state.
%
% Prints model dimensions, the objective function, and — when a flux vector
% is supplied — the objective value plus all net exchange fluxes partitioned
% into uptake and secretion. Inspired by COBRApy's model.summary().
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure.
%
% Name-Value Arguments
% --------------------
% fluxes : double
%     flux vector (length = numel(model.rxns)). When provided, shows the
%     objective value and the exchange-flux summary (default []).
% cutoff : double
%     minimum absolute net exchange flux to include (default 1e-8).
% nMax : double
%     maximum exchange reactions displayed per direction (default 20).
%
% Examples
% --------
%     modelSummary(model)
%     modelSummary(model, 'fluxes', sol.x)
%     modelSummary(model, 'fluxes', sol.x, 'nMax', 10)

p = parseRAVENargs(varargin, {'fluxes',[]; 'cutoff',1e-8; 'nMax',20});
fluxes = p.fluxes;
cutoff = p.cutoff;
nMax   = p.nMax;

nRxns  = numel(model.rxns);
nMets  = numel(model.mets);
nGenes = 0; if isfield(model,'genes'), nGenes = numel(model.genes); end
nComps = 0; if isfield(model,'comps'), nComps = numel(model.comps); end

W   = 62;
bar = repmat('=', 1, W);
sep = repmat('-', 1, W);

fprintf('%s\n', bar);
modelID = ''; if isfield(model,'id') && ~isempty(model.id), modelID = model.id; end
if isfield(model,'name') && ~isempty(model.name)
    fprintf('  %s  (%s)\n', model.name, modelID);
elseif ~isempty(modelID)
    fprintf('  %s\n', modelID);
else
    fprintf('  (unnamed model)\n');
end
fprintf('%s\n', bar);
fprintf('  %-15s %5d    %-15s %5d\n', 'Reactions:',  nRxns, 'Metabolites:',  nMets);
fprintf('  %-15s %5d    %-15s %5d\n', 'Genes:',      nGenes,'Compartments:', nComps);
fprintf('%s\n', sep);

% Objective
objIdx = find(model.c ~= 0);
if isempty(objIdx)
    fprintf('  Objective: (none)\n');
else
    fprintf('  Objective:\n');
    for i = 1:numel(objIdx)
        j = objIdx(i);
        rName = '';
        if isfield(model,'rxnNames') && j <= numel(model.rxnNames) && ~isempty(model.rxnNames{j})
            rName = [' (' model.rxnNames{j} ')'];
        end
        fprintf('    %+.4g \xD7 %s%s\n', model.c(j), model.rxns{j}, rName);
    end
end

if isempty(fluxes)
    fprintf('%s\n', bar);
    return
end

objVal = model.c' * fluxes;
fprintf('  Objective value: %.6g\n', objVal);
fprintf('%s\n', sep);

% Classify exchange fluxes
[~, exchIdx] = getExchangeRxns(model);

upRows  = {};   % {metLabel, rxnID, absFlux}
outRows = {};

for i = 1:numel(exchIdx)
    f = fluxes(exchIdx(i));
    if abs(f) < cutoff, continue; end
    col = model.S(:, exchIdx(i));
    m = find(col, 1);           % exchange reactions have exactly one metabolite
    if isempty(m), continue; end
    netProd = full(col(m)) * f; % positive = metabolite secreted; full() avoids sparse scalar in cell
    mLabel  = metLabel_(model, m);
    rID     = model.rxns{exchIdx(i)};
    if netProd > 0
        outRows{end+1} = {mLabel, rID, netProd};   %#ok<AGROW>
    else
        upRows{end+1}  = {mLabel, rID, -netProd};  %#ok<AGROW>
    end
end

printSection_('  UPTAKE',    upRows,  nMax, W);
printSection_('  SECRETION', outRows, nMax, W);
fprintf('%s\n', bar);
end

%--------------------------------------------------------------------------
function printSection_(title, rows, nMax, W)
    fprintf('\n%s\n', title);
    fprintf('  %-26s  %-20s  %10s\n', 'Metabolite', 'Reaction', 'Flux');
    fprintf('  %s\n', repmat('-', 1, W-2));
    if isempty(rows)
        fprintf('  (none)\n');
        return
    end
    vals = cellfun(@(r) r{3}, rows);
    [~, si] = sort(vals, 'descend');
    rows = rows(si);
    for i = 1:min(numel(rows), nMax)
        fprintf('  %-26s  %-20s  %10.4g\n', rows{i}{1}, rows{i}{2}, rows{i}{3});
    end
    if numel(rows) > nMax
        fprintf('  ... (%d more, use nMax to show all)\n', numel(rows) - nMax);
    end
end

function s = metLabel_(model, metIdx)
    if isfield(model,'metNames') && ~isempty(model.metNames{metIdx})
        s = model.metNames{metIdx};
    else
        s = model.mets{metIdx};
    end
    if numel(s) > 26, s = [s(1:23) '...']; end
end
