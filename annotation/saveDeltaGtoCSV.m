function saveDeltaGtoCSV(model, varargin)
% saveDeltaGtoCSV  Persist metDeltaG and rxnDeltaG to CSV files.
%
% Persist model.metDeltaG and model.rxnDeltaG to project CSV files.
% Counterpart of loadDeltaGfromCSV and the upstream version of yeast-GEM's
% saveDeltaG. Mirrors raven_python.annotation.save_delta_g_csv.
%
% Each CSV gets two columns: identifier, deltaG. Rows are written in model
% order (one row per entity); identifiers without a matching field get
% NaN. Pass an empty string for metCsv or rxnCsv to skip that side.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% metCsv : char, optional
%     Output path for the metabolite ΔG CSV, or '' to skip.
% rxnCsv : char, optional
%     Output path for the reaction ΔG CSV, or '' to skip.
% verbose : logical, optional
%     Print "wrote ..." per file (default false).
%
% Examples
% --------
%     saveDeltaGtoCSV(model, ...
%         'data/databases/model_metDeltaG.csv', ...
%         'data/databases/model_rxnDeltaG.csv');

p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''; 'verbose',false});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;
verbose=p.verbose;

if ~isempty(metCsv)
    if ~isfield(model, 'metDeltaG')
        if verbose
            fprintf('No metDeltaG field found, %s will not be changed.\n', metCsv);
        end
    else
        metG = array2table([model.mets, num2cell(model.metDeltaG)]);
        writetable(metG, metCsv);
        if verbose
            fprintf('Wrote %s\n', metCsv);
        end
    end
end

if ~isempty(rxnCsv)
    if ~isfield(model, 'rxnDeltaG')
        if verbose
            fprintf('No rxnDeltaG field found, %s will not be changed.\n', rxnCsv);
        end
    else
        rxnG = array2table([model.rxns, num2cell(model.rxnDeltaG)]);
        writetable(rxnG, rxnCsv);
        if verbose
            fprintf('Wrote %s\n', rxnCsv);
        end
    end
end
end
