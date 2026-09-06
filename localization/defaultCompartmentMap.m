function map = defaultCompartmentMap()
% defaultCompartmentMap  Default predictor/database compartment label -> model compartment id.
%
% Returns a containers.Map from lower-case predictor / database compartment labels (as used by
% DeepLoc, MULocDeep, COMPARTMENTS and UniProt) to compartment ids, tuned for yeast/fungal
% models (e.g. yeast-GEM codes). Pass it (or your own map) to parseScores / getUniProtScores as
% 'compartmentMap'.
%
% The table below is also the worked example of how such a map is built, so read it as three
% deliberate decisions rather than a lookup list:
%
% - **Which compartments exist.** A label absent from the map is dropped along with its score
%   column, so listing only the ids your model has is how unwanted compartments are excluded.
%   'plastid' is absent because fungi lack one; a plant model would add it.
% - **Which labels are the same thing.** Predictors disagree on wording, so 'cytoplasm' and
%   'cytosol' both map to 'c', and 'mitochondrion' / 'mitochondria' / 'mitochondrial' all map
%   to 'm'. Merged columns are combined by maximum, not sum.
% - **Which distinct compartments to collapse.** 'lysosome' maps to 'v' not because it is a
%   vacuole but because a yeast model has no lysosome and the vacuole is its functional
%   equivalent. Change this if your model separates them.
%
% Extend a copy rather than editing this function, e.g.
%
%     map = defaultCompartmentMap;
%     map('cell wall') = 'ce';
%     GSS = parseScores(file, 'compartmentMap', map);
%
% Returns
% -------
% map : containers.Map
%     label (char, lower case) -> compartment id (char).
%
% See also
% --------
% parseScores, getUniProtScores, predictLocalization

labels = {'cytoplasm','cytosol','nucleus','nucleoplasm','mitochondrion','mitochondria', ...
          'mitochondrial','peroxisome','endoplasmic reticulum','golgi apparatus','golgi', ...
          'vacuole','lysosome/vacuole','lysosome','extracellular','extracellular space', ...
          'extracellular region','secreted','cell membrane','plasma membrane','cell envelope', ...
          'lipid particle','lipid droplet'};
ids    = {'c','c','n','n','m','m', ...
          'm','p','er','g','g', ...
          'v','v','v','e','e', ...
          'e','e','ce','ce','ce', ...
          'lp','lp'};
map = containers.Map(labels, ids);
end
