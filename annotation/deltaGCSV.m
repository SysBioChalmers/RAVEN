function model = deltaGCSV(model, direction, varargin)
% deltaGCSV  Load or save metDeltaG and rxnDeltaG via CSV files.
%
% Consolidates loadDeltaGfromCSV and saveDeltaGtoCSV into a single function.
% Mirrors raven_python.annotation.load_delta_g_csv /save_delta_g_csv.
%
% Each CSV is a two-column table: identifier, deltaG. On load, rows whose
% identifier is not in the model are silently skipped. Pass '' for either
% metCsv or rxnCsv to skip that side.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% direction : char
%     'load' to populate model fields from CSV; 'save' to write them out.
%
% Name-Value Arguments
% --------------------
% metCsv : char
%     Path to the metabolite deltaG CSV, or '' to skip.
% rxnCsv : char
%     Path to the reaction deltaG CSV, or '' to skip.
% verbose : logical
%     Print a confirmation line per file written (save only; default false).
%
% Returns
% -------
% model : struct
%     Updated model struct (unchanged for save direction).
%
% Examples
% --------
%     model = deltaGCSV(model, 'load', 'data/met_dG.csv', 'data/rxn_dG.csv');
%     deltaGCSV(model, 'save', 'data/met_dG.csv', 'data/rxn_dG.csv');

direction=char(direction);
p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''; 'verbose',false});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;
verbose=p.verbose;

switch direction
    case 'load'
        model=loadField(model, 'metDeltaG', model.mets, metCsv);
        model=loadField(model, 'rxnDeltaG', model.rxns, rxnCsv);
    case 'save'
        saveField(model, 'metDeltaG', model.mets, metCsv, verbose);
        saveField(model, 'rxnDeltaG', model.rxns, rxnCsv, verbose);
    otherwise
        error('direction must be ''load'' or ''save''');
end
end

function model=loadField(model, fieldname, ids, csvPath)
if isempty(csvPath); return; end
if isfield(model, fieldname)
    disp(['Existing ' fieldname ' field will be overwritten.'])
else
    model.(fieldname)=nan(numel(ids),1);
end
G=readtable(csvPath);
[a, b]=ismember(ids, G.(G.Properties.VariableNames{1}));
model.(fieldname)(a)=G.(G.Properties.VariableNames{2})(b(a));
if any(~a)
    fprintf(['Not all identifiers are matched to %s; the latter\n' ...
             'file might have to be supplemented with deltaG values for new entries.\n'], ...
             csvPath);
end
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
