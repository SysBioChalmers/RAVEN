function model = scaleBiomassFraction(model, biomassConfig, componentName, newValue, balanceOut)
% scaleBiomassFraction  Rescale a biomass component to a target value.
%
% Rescale a biomass component to a target g/gDW value, optionally
% balancing a second component so the total biomass mass stays at 1 g/gDW.
% Mirrors raven_python.biomass.scale_biomass and yeast-GEM's legacy
% scaleBioMass.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassConfig : struct
%     Struct (see getBiomassFractions).
% componentName : char
%     Component to rescale.
% newValue : double
%     Target fraction in g/gDW.
% balanceOut : char, optional
%     Second component name to adjust so the biomass total remains 1
%     g/gDW. Empty / omit to skip balancing.
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = scaleBiomassFraction(model, biomassConfig, 'protein', 0.5, 'carbohydrate');
%
% See also
% --------
% getBiomassFractions

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
