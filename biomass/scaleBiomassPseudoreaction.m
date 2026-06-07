function model = scaleBiomassPseudoreaction(model, biomassConfig, componentName, factor)
% scaleBiomassPseudoreaction
%   Multiply the substrate coefficients of one biomass component
%   pseudoreaction by `factor` and rebalance H+ to preserve charge
%   neutrality. Mirrors raven_python.biomass.rescale_pseudoreaction
%   and yeast-GEM's legacy rescalePseudoReaction.
%
%   "Substrate" means every metabolite in the pseudoreaction whose
%   metabolite name does NOT match the component name (the
%   component's product is left untouched). After rescaling, the
%   coefficient of biomassConfig.proton_met is recomputed so the
%   pseudoreaction's total ionic charge sums to zero.
%
%   Inputs:
%       model           RAVEN model struct.
%       biomassConfig   struct (see getBiomassFractions).
%       componentName   Name of the component to rescale (must match
%                       biomassConfig.components{i}.name for some i,
%                       AND be the model.metNames of the produced
%                       metabolite in the matching pseudoreaction).
%       factor          Multiplicative factor.
%
%   Output:
%       model           Modified model.
%
% Usage: model = scaleBiomassPseudoreaction(model, biomassConfig, 'protein', 0.9)

comp = findComponent(biomassConfig, componentName);
rxnPos = find(strcmp(model.rxnNames, comp.pseudoreaction_name));
if isempty(rxnPos)
    error('scaleBiomassPseudoreaction:missingPseudoreaction', ...
        'No reaction named %s in model.', comp.pseudoreaction_name);
end

for i = 1:length(model.mets)
    S_ir = model.S(i, rxnPos);
    isProd = strcmp(model.metNames{i}, componentName);
    if S_ir ~= 0 && ~isProd
        model.S(i, rxnPos) = factor * S_ir;
    end
end

% Rebalance H+ to keep charge neutrality.
Hc = find(strcmp(model.mets, biomassConfig.proton_met));
model.S(Hc, rxnPos) = 0;
model.S(Hc, rxnPos) = -sum(model.S(:, rxnPos) .* model.metCharges, 'omitnan');
end

function comp = findComponent(cfg, name)
for i = 1:numel(cfg.components)
    if strcmp(cfg.components{i}.name, name)
        comp = cfg.components{i};
        return;
    end
end
error('scaleBiomassPseudoreaction:unknownComponent', ...
    'biomassConfig has no component named %s', name);
end
