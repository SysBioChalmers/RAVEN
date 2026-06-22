function map = defaultCompartmentMap()
% defaultCompartmentMap  Default predictor/database compartment label -> model compartment id.
%
% Returns a containers.Map from lower-case predictor / database compartment labels (as used by
% DeepLoc, MULocDeep, COMPARTMENTS and UniProt) to compartment ids, tuned for yeast/fungal
% models (e.g. yeast-GEM codes). Pass it (or your own map) to parseScores / getUniProtScores as
% 'compartmentMap' to rename and merge compartments; labels not listed here (e.g. plastid, which
% fungi lack) are dropped.
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
