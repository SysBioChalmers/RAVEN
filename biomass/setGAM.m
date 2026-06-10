function model = setGAM(model, value, biomassRxn, cofactorMetNames, ngamRxn, ngamValue)
% setGAM  Set the growth-associated maintenance (GAM) coefficient.
%
% Set the growth-associated maintenance (GAM) coefficient in the biomass
% pseudoreaction, and optionally fix the non-growth maintenance (NGAM)
% reaction's bounds. Mirrors raven_python.biomass.set_gam and yeast-GEM's
% legacy changeGAM.
%
% For every metabolite in the biomass pseudoreaction whose model.metNames
% entry is in cofactorMetNames, the stoichiometric coefficient is set to
% ±value preserving the sign of the current coefficient. Yeast-GEM scales
% ATP, ADP, H2O, H+ and phosphate (with ATP and H2O on the substrate side,
% ADP / H+ / phosphate on the product side).
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% value : double
%     New GAM value (mmol ATP / gDW per growth unit).
% biomassRxn : char
%     Reaction id of the biomass pseudoreaction.
% cofactorMetNames : cell
%     Cell array of metabolite NAMES (not IDs) to rescale, e.g.
%     {'ATP','ADP','H2O','H+','phosphate'}.
% ngamRxn : char, optional
%     NGAM reaction id. Required when ngamValue is supplied.
% ngamValue : double, optional
%     NGAM flux to fix. Sets the NGAM reaction's bounds to (ngamValue,
%     ngamValue).
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = setGAM(model, 80, 'r_4041', {'ATP','ADP','H2O','H+','phosphate'});

if nargin < 4
    error('setGAM:missingArgs', ...
        'biomassRxn and cofactorMetNames are required.');
end

bioPos = strcmp(model.rxns, biomassRxn);
if ~any(bioPos)
    error('setGAM:missingBiomassRxn', ...
        'Reaction %s not found in model.', biomassRxn);
end

for i = 1:length(model.mets)
    S_ix = model.S(i, bioPos);
    isCofactor = any(strcmp(cofactorMetNames, model.metNames{i}));
    if S_ix ~= 0 && isCofactor
        model.S(i, bioPos) = sign(S_ix) * value;
    end
end

if nargin >= 6 && ~isempty(ngamRxn) && ~isempty(ngamValue)
    model = setParam(model, 'eq', ngamRxn, ngamValue);
end
end
