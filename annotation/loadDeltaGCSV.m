function model = loadDeltaGCSV(model, varargin)
% loadDeltaGCSV  Load metDeltaG and rxnDeltaG from CSV files.
%
% Mirrors raven_toolbox.annotation.load_delta_g_csv.
%
% Each CSV is a two-column table: identifier, deltaG. On load, rows whose
% identifier is not in the model are silently skipped, and every matched
% value is stored exactly as it appears in the CSV -- including yeast-GEM's
% own "no measurement" placeholder (10000000.0). Callers that need to treat
% a particular value as absent should filter it themselves after loading;
% this function does not interpret CSV values. Pass '' for either metCsv
% or rxnCsv to skip that side.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
%
% Name-Value Arguments
% --------------------
% metCsv : char
%     Path to the metabolite deltaG CSV, or '' to skip.
% rxnCsv : char
%     Path to the reaction deltaG CSV, or '' to skip.
%
% Returns
% -------
% model : struct
%     Updated model struct.
%
% Examples
% --------
%     model = loadDeltaGCSV(model, 'data/met_dG.csv', 'data/rxn_dG.csv');
%     model = loadDeltaGCSV(model, 'metCsv', 'data/met_dG.csv');

p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;

model=loadField(model, 'metDeltaG', model.mets, metCsv);
model=loadField(model, 'rxnDeltaG', model.rxns, rxnCsv);
end

function model=loadField(model, fieldname, ids, csvPath)
if isempty(csvPath); return; end
if isfield(model, fieldname)
    disp(['Existing ' fieldname ' field will be overwritten.'])
else
    model.(fieldname)=nan(numel(ids),1);
end
G=readtable(csvPath);
values=G.(G.Properties.VariableNames{2});
[a, b]=ismember(ids, G.(G.Properties.VariableNames{1}));
model.(fieldname)(a)=values(b(a));
if any(~a)
    fprintf(['Not all identifiers are matched to %s; the latter\n' ...
             'file might have to be supplemented with deltaG values for new entries.\n'], ...
             csvPath);
end
end
