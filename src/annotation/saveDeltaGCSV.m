function saveDeltaGCSV(model, varargin)
% saveDeltaGCSV  Save metDeltaG and rxnDeltaG to CSV files.
%
% Writes a two-column table (identifier, deltaG) per side, one row per
% entity, in model order. An entity with no value stored is written as
% NaN, unchanged from whatever the model field holds -- this function does
% not interpret or substitute values. Pass '' for either metCsv or rxnCsv
% to skip that side.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
%
% Name-Value Arguments
% --------------------
% metCsv : char
%     Path to write the metabolite deltaG CSV to, or '' to skip.
% rxnCsv : char
%     Path to write the reaction deltaG CSV to, or '' to skip.
% verbose : logical
%     Print a confirmation line per file written (default false).
%
% Examples
% --------
%     saveDeltaGCSV(model, 'data/met_dG.csv', 'data/rxn_dG.csv');
%     saveDeltaGCSV(model, 'metCsv', 'data/met_dG.csv', 'verbose', true);

p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''; 'verbose',false});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;
verbose=p.verbose;

saveField(model, 'metDeltaG', model.mets, metCsv, verbose);
saveField(model, 'rxnDeltaG', model.rxns, rxnCsv, verbose);
end

function saveField(model, fieldname, ids, csvPath, verbose)
if isempty(csvPath); return; end
if ~isfield(model, fieldname)
    if verbose
        fprintf('No %s field found, %s will not be changed.\n', fieldname, csvPath);
    end
else
    G=array2table([ids, num2cell(model.(fieldname))]);
    writetable(G, csvPath);
    if verbose
        fprintf('Wrote %s\n', csvPath);
    end
end
end
