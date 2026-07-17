function report = diffModels(modelA, modelB, varargin)
% diffModels  Report the semantic differences between two models.
%
% Counterpart of raven_toolbox.comparison.diff.diff_models. Compares two
% models entry by entry, keyed by identifier, and returns a report listing
% every difference found. Reactions, metabolites and genes are matched by
% id; for the ids present in both models the reaction stoichiometry, bounds,
% objective coefficient, grRule, EC codes, and the metabolite formula,
% charge and compartment are compared. This is the structured
% added/removed/changed diff that compareMultipleModels /
% compareRxnsGenesMetsComps (which report overlap counts and distances) do
% not provide.
%
% grRules are compared for logical equality rather than as strings: each
% rule is expanded to disjunctive normal form and the genes within each
% isozyme, and the isozymes themselves, are sorted before comparison, so
% "a and b" equals "b and a" and "(a or b)" equals "(b or a)". This is
% stronger than a case/whitespace string match.
%
% Parameters
% ----------
% modelA : struct
%     the first model.
% modelB : struct
%     the second model.
%
% Name-Value Arguments
% --------------------
% stoichTol : double
%     absolute tolerance for stoichiometric coefficient differences
%     (default 1e-9).
% maxPerCategory : double
%     the maximum number of differences reported per category before the
%     rest are truncated with a single summary line (default 50).
% printResults : logical
%     print the report to the command window (default false).
%
% Returns
% -------
% report : struct
%     difference report with fields:
%
%     - equal : true if no differences were found.
%     - differences : cell array of human-readable difference strings, each
%       naming the category and the entity affected.
%     - rxnsOnlyInA / rxnsOnlyInB : reaction ids present in only one model.
%     - metsOnlyInA / metsOnlyInB : metabolite ids present in only one model.
%     - genesOnlyInA / genesOnlyInB : gene ids present in only one model.
%
% Usage: report = diffModels(modelA, modelB, ...)

p = parseRAVENargs(varargin, {'stoichTol', 1e-9; 'maxPerCategory', 50; ...
    'printResults', false});
stoichTol = p.stoichTol;
maxPerCategory = p.maxPerCategory;
printResults = p.printResults;

diffs = {};                 % collected difference strings
counts = struct();          % per-category counters, for truncation

report = struct();
report.equal = true;
report.differences = {};

% --- Identifier sets --------------------------------------------------
[report.rxnsOnlyInA, report.rxnsOnlyInB, commonRxns] = idSets(modelA.rxns, modelB.rxns);
[report.metsOnlyInA, report.metsOnlyInB, commonMets] = idSets(modelA.mets, modelB.mets);
if isfield(modelA,'genes') && isfield(modelB,'genes')
    [report.genesOnlyInA, report.genesOnlyInB] = idSets(modelA.genes, modelB.genes);
else
    report.genesOnlyInA = {}; report.genesOnlyInB = {};
end

reportSet('reactions', report.rxnsOnlyInA, report.rxnsOnlyInB);
reportSet('metabolites', report.metsOnlyInA, report.metsOnlyInB);
reportSet('genes', report.genesOnlyInA, report.genesOnlyInB);

% --- Reaction content, on the shared ids ------------------------------
[~, ia] = ismember(commonRxns, modelA.rxns);
[~, ib] = ismember(commonRxns, modelB.rxns);
for k = 1:numel(commonRxns)
    rxn = commonRxns{k};
    ja = ia(k); jb = ib(k);

    % Stoichiometry: compare the metabolite id -> coefficient maps.
    [ma, ca] = rxnStoich(modelA, ja);
    [mb, cb] = rxnStoich(modelB, jb);
    if ~isequal(sort(ma), sort(mb))
        push('stoich', sprintf('%s: stoichiometry involves different metabolites', rxn));
    else
        [~, locb] = ismember(ma, mb);
        for m = 1:numel(ma)
            if abs(ca(m) - cb(locb(m))) > stoichTol
                push('stoich', sprintf('%s: coef[%s] A=%g B=%g', rxn, ma{m}, ca(m), cb(locb(m))));
            end
        end
    end

    % Bounds.
    if modelA.lb(ja) ~= modelB.lb(jb) || modelA.ub(ja) ~= modelB.ub(jb)
        push('bounds', sprintf('%s: bounds A=[%g,%g] B=[%g,%g]', rxn, ...
            modelA.lb(ja), modelA.ub(ja), modelB.lb(jb), modelB.ub(jb)));
    end

    % Objective coefficient.
    if fieldAt(modelA,'c',ja,0) ~= fieldAt(modelB,'c',jb,0)
        push('objective', sprintf('%s: objective A=%g B=%g', rxn, ...
            fieldAt(modelA,'c',ja,0), fieldAt(modelB,'c',jb,0)));
    end

    % grRule, compared as logic rather than text.
    ga = canonicalGpr(cellAt(modelA,'grRules',ja));
    gb = canonicalGpr(cellAt(modelB,'grRules',jb));
    if ~strcmp(ga, gb)
        push('gpr', sprintf('%s: grRule A="%s" B="%s"', rxn, ...
            cellAt(modelA,'grRules',ja), cellAt(modelB,'grRules',jb)));
    end

    % EC codes (RAVEN's first-class annotation field).
    ea = cellAt(modelA,'eccodes',ja);
    eb = cellAt(modelB,'eccodes',jb);
    if ~strcmp(ea, eb)
        push('eccodes', sprintf('%s: eccodes A="%s" B="%s"', rxn, ea, eb));
    end
end

% --- Metabolite content, on the shared ids ----------------------------
[~, ia] = ismember(commonMets, modelA.mets);
[~, ib] = ismember(commonMets, modelB.mets);
for k = 1:numel(commonMets)
    met = commonMets{k};
    ja = ia(k); jb = ib(k);

    fa = cellAt(modelA,'metFormulas',ja);
    fb = cellAt(modelB,'metFormulas',jb);
    if ~strcmp(fa, fb)
        push('formula', sprintf('%s: formula A="%s" B="%s"', met, fa, fb));
    end

    if ~isequaln(fieldAt(modelA,'metCharges',ja,NaN), fieldAt(modelB,'metCharges',jb,NaN))
        push('charge', sprintf('%s: charge A=%g B=%g', met, ...
            fieldAt(modelA,'metCharges',ja,NaN), fieldAt(modelB,'metCharges',jb,NaN)));
    end

    cpA = metComp(modelA, ja);
    cpB = metComp(modelB, jb);
    if ~strcmp(cpA, cpB)
        push('compartment', sprintf('%s: compartment A="%s" B="%s"', met, cpA, cpB));
    end
end

report.equal = isempty(diffs);
report.differences = diffs(:);

if printResults
    if report.equal
        fprintf('Models are semantically equal.\n');
    else
        fprintf('Models differ (%d differences):\n', numel(diffs));
        fprintf('  - %s\n', diffs{:});
    end
end

    % --- nested helpers (share diffs / counts / maxPerCategory) --------
    function push(category, message)
        if ~isfield(counts, category)
            counts.(category) = 0;
        end
        counts.(category) = counts.(category) + 1;
        if counts.(category) <= maxPerCategory
            diffs{end+1,1} = message; %#ok<AGROW>
        elseif counts.(category) == maxPerCategory + 1
            diffs{end+1,1} = sprintf('... (%s category truncated at %d entries)', ...
                category, maxPerCategory); %#ok<AGROW>
        end
    end

    function reportSet(label, onlyA, onlyB)
        if ~isempty(onlyA)
            push([label '_only_in_a'], sprintf('%d %s only in A: %s', ...
                numel(onlyA), label, joinTruncated(onlyA)));
        end
        if ~isempty(onlyB)
            push([label '_only_in_b'], sprintf('%d %s only in B: %s', ...
                numel(onlyB), label, joinTruncated(onlyB)));
        end
    end
end

% --- local helpers ----------------------------------------------------
function [onlyA, onlyB, common] = idSets(a, b)
onlyA = setdiff(a, b);
onlyB = setdiff(b, a);
common = sort(intersect(a, b));
onlyA = sort(onlyA(:))'; onlyB = sort(onlyB(:))';
onlyA = onlyA(:); onlyB = onlyB(:);
end

function [mets, coeffs] = rxnStoich(model, j)
nz = find(model.S(:, j));
mets = model.mets(nz);
mets = mets(:)';
coeffs = full(model.S(nz, j))';
end

function v = fieldAt(model, field, j, default)
if isfield(model, field) && numel(model.(field)) >= j
    v = model.(field)(j);
    if iscell(v), v = v{1}; end
    if isempty(v), v = default; end
else
    v = default;
end
end

function s = cellAt(model, field, j)
% A char scalar from a cell-array field, '' when absent or empty.
if isfield(model, field) && numel(model.(field)) >= j
    s = model.(field){j};
    if isempty(s), s = ''; end
else
    s = '';
end
end

function cp = metComp(model, j)
if isfield(model, 'metComps') && isfield(model, 'comps')
    cp = model.comps{model.metComps(j)};
else
    cp = '';
end
end

function s = canonicalGpr(rule)
% A canonical string for a grRule under which isozymes and their subunits
% are order-insensitive. Empty rule -> ''.
if isempty(rule)
    s = '';
    return;
end
clauses = grRuleToDNF(rule);
parts = cell(1, numel(clauses));
for i = 1:numel(clauses)
    genes = sort(unique(clauses{i}));
    parts{i} = strjoin(genes(:)', ' and ');
end
parts = sort(parts);
s = strjoin(parts, ' or ');
end

function s = joinTruncated(ids)
% Join up to the first 10 ids, appending a count of the remainder.
ids = ids(:)';
if numel(ids) > 10
    s = [strjoin(ids(1:10), ', ') sprintf(', ... (%d more)', numel(ids) - 10)];
else
    s = strjoin(ids, ', ');
end
end
