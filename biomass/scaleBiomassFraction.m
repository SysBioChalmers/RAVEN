function model = scaleBiomassFraction(model, biomassConfig, componentName, newValue, balanceOut)
% scaleBiomassFraction
%   Rescale a biomass component to a target g/gDW value, optionally
%   balancing a second component so the total biomass mass stays at
%   1 g/gDW. Mirrors raven_python.biomass.scale_biomass and yeast-GEM's
%   legacy scaleBioMass.
%
%   Inputs:
%       model           RAVEN model struct.
%       biomassConfig   struct (see getBiomassFractions).
%       componentName   Component to rescale.
%       newValue        Target fraction in g/gDW.
%       balanceOut      (opt) Second component name to adjust so the
%                       biomass total remains 1 g/gDW. Empty / omit
%                       to skip balancing.
%
%   Output:
%       model           Modified model.
%
% Usage: model = scaleBiomassFraction(model, biomassConfig, 'protein', 0.5, 'carbohydrate')

if nargin < 5
    balanceOut = '';
end

fractions = getBiomassFractions(model, biomassConfig);
if ~isfield(fractions, componentName)
    error('scaleBiomassFraction:unknownComponent', ...
        'biomassConfig has no component named %s', componentName);
end
current = fractions.(componentName);
if current == 0
    error('scaleBiomassFraction:zeroCurrent', ...
        ['Cannot scale %s to %g: current fraction is 0 ' ...
         '(pseudoreaction missing or empty).'], componentName, newValue);
end
factor = newValue / current;
model = scaleBiomassPseudoreaction(model, biomassConfig, componentName, factor);

if ~isempty(balanceOut)
    fractions = getBiomassFractions(model, biomassConfig);
    total = fractions.total;
    balanceCurrent = fractions.(balanceOut);
    if balanceCurrent == 0
        return;
    end
    balanceFactor = (balanceCurrent + (1 - total)) / balanceCurrent;
    model = scaleBiomassPseudoreaction(model, biomassConfig, balanceOut, balanceFactor);
end
end
