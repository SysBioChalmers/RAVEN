function model = loadDeltaGfromCSV(model, varargin)
% loadDeltaGfromCSV  Populate metDeltaG and rxnDeltaG from CSV files.
%
% Populate model.metDeltaG and model.rxnDeltaG from project CSV files.
% Mirrors raven_python.annotation.load_delta_g_csv and is the upstream
% version of yeast-GEM's loadDeltaG.
%
% Each CSV is a two-column table: identifier, deltaG. Rows whose
% identifier doesn't appear in the model are silently skipped. Pass an
% empty string ('') for either argument to skip that side.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
%
% Name-Value Arguments
% --------------------
% metCsv : char
%     Path to metabolite ΔG CSV (id, ΔG), or '' to skip.
% rxnCsv : char
%     Path to reaction ΔG CSV (id, ΔG), or '' to skip.
%
% Returns
% -------
% model : struct
%     Model with metDeltaG and/or rxnDeltaG fields added.
%
% Examples
% --------
%     model = loadDeltaGfromCSV(model, ...
%         'data/databases/model_metDeltaG.csv', ...
%         'data/databases/model_rxnDeltaG.csv');

p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;

if ~isempty(metCsv)
    if isfield(model, 'metDeltaG')
        disp('Existing metDeltaG field will be overwritten.')
    else
        model.metDeltaG = nan(numel(model.mets), 1);
    end
    metG = readtable(metCsv);
    [a, b] = ismember(model.mets, metG.(metG.Properties.VariableNames{1}));
    model.metDeltaG(a) = metG.(metG.Properties.VariableNames{2})(b(a));
    if any(~a)
        fprintf(['Not all metabolite identifiers are matched to %s; the latter\n' ...
                 'file might have to be supplemented with deltaG values for new metabolites.\n'], ...
                 metCsv);
    end
end

if ~isempty(rxnCsv)
    if isfield(model, 'rxnDeltaG')
        disp('Existing rxnDeltaG field will be overwritten.')
    else
        model.rxnDeltaG = nan(numel(model.rxns), 1);
    end
    rxnG = readtable(rxnCsv);
    [a, b] = ismember(model.rxns, rxnG.(rxnG.Properties.VariableNames{1}));
    model.rxnDeltaG(a) = rxnG.(rxnG.Properties.VariableNames{2})(b(a));
    if any(~a)
        fprintf(['Not all reaction identifiers are matched to %s; the latter\n' ...
                 'file might have to be supplemented with deltaG values for new reactions.\n'], ...
                 rxnCsv);
    end
end
end
