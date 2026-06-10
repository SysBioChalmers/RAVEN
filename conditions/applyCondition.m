function model = applyCondition(model, condition)
% applyCondition  Apply a deterministic condition to a model.
%
% Apply a deterministic "condition" to a model: a prelude that resets
% exchange bounds, optional metabolite removals + automatic charge
% rebalancing of a pseudoreaction, optional biomass-stoichiometry delta,
% and a per-reaction bounds diff. The schema is intentionally narrow so a
% condition can be reviewed as data.
%
% Yeast-GEM was the first consumer; the same schema works for any GEM that
% keeps its condition presets as data rather than as code.
% Project-specific extensions (e.g. yeast-GEM's amino_acid_ratio step that
% rewrites a protein pseudoreaction's stoichiometry from a side-car TSV)
% are handled by the caller before / after this function — kept
% upstream-narrow on purpose.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% condition : char or struct
%     Either a path to a YAML condition file or a struct already produced
%     by parseYAML. The expected schema (all keys optional):
%
%         prelude:
%           reset_exchanges: out      % truthy -> reset all
%
%         cofactor_pseudoreaction:
%           rxn_id: r_4598
%           remove_mets:
%             - { met: s_3714 }
%           charge_balance_met: s_0794
%
%         biomass_stoichiometry_delta:
%           rxn_id: r_4041
%           add:
%             - { met: s_0689, coef:  0.08 }
%             - { met: s_0687, coef: -0.08 }
%             - { met: s_0794, coef: -0.16 }
%
%         bounds:
%           - { rxn: r_1654, lb: -1000 }
%           - { rxn: r_1992, lb: 0 }
%           - { rxn: r_1663, lb: 0, ub: 0 }
%
%         expected_uptake_count: 15
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = applyCondition(model, 'data/conditions/anaerobic.yml');
%     model = applyCondition(model, parseYAML('data/conditions/anaerobic.yml'));
%
% See also
% --------
% parseYAML

if ischar(condition) || isstring(condition)
    cond = parseYAML(char(condition));
elseif isstruct(condition)
    cond = condition;
else
    error('applyCondition:invalidCondition', ...
        'condition must be a YAML file path or a struct.');
end

% --- Step 1: prelude ---------------------------------------------------
if isfield(cond, 'prelude') && isfield(cond.prelude, 'reset_exchanges')
    [~, exchangeRxns] = getExchangeRxns(model, cond.prelude.reset_exchanges);
    model.lb(exchangeRxns) = 0;
    model.ub(exchangeRxns) = 1000;
end

% --- Step 2: cofactor pseudoreaction edits ----------------------------
if isfield(cond, 'cofactor_pseudoreaction')
    cp = cond.cofactor_pseudoreaction;
    cofacIdx = getIndexes(model, cp.rxn_id, 'rxns');
    if isfield(cp, 'remove_mets')
        for i = 1:numel(cp.remove_mets)
            metIdx = getIndexes(model, cp.remove_mets{i}.met, 'mets');
            model.S(metIdx, cofacIdx) = 0;
        end
    end
    if isfield(cp, 'charge_balance_met')
        balanceIdx = find(strcmp(model.mets, cp.charge_balance_met));
        model.S(balanceIdx, cofacIdx) = 0;
        model.S(balanceIdx, cofacIdx) = ...
            -sum(model.S(:, cofacIdx) .* model.metCharges, 'omitnan');
    end
end

% --- Step 3: biomass stoichiometry delta ------------------------------
if isfield(cond, 'biomass_stoichiometry_delta')
    delta = cond.biomass_stoichiometry_delta;
    bioIdx = getIndexes(model, delta.rxn_id, 'rxns');
    if isfield(delta, 'add')
        for i = 1:numel(delta.add)
            entry = delta.add{i};
            metIdx = getIndexes(model, entry.met, 'mets');
            model.S(metIdx, bioIdx) = full(model.S(metIdx, bioIdx)) + entry.coef;
        end
    end
end

% --- Step 4: bounds ---------------------------------------------------
nUptake = 0;
if isfield(cond, 'bounds')
    for i = 1:numel(cond.bounds)
        b = cond.bounds{i};
        rxnIdx = find(strcmp(model.rxns, b.rxn));
        if isempty(rxnIdx)
            warning('applyCondition:missingRxn', ...
                'Reaction %s not found in model; skipping.', b.rxn);
            continue;
        end
        if isfield(b, 'lb')
            model.lb(rxnIdx) = b.lb;
            if b.lb == -1000
                nUptake = nUptake + 1;
            end
        end
        if isfield(b, 'ub')
            model.ub(rxnIdx) = b.ub;
        end
    end
end

% --- Step 5: uptake sanity check --------------------------------------
if isfield(cond, 'expected_uptake_count')
    if nUptake ~= cond.expected_uptake_count
        warning('applyCondition:uptakeCountMismatch', ...
            'Expected %d uptake reactions, applied %d. Some may be missing from the model.', ...
            cond.expected_uptake_count, nUptake);
    end
end

end
