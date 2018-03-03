function model=replaceMets(model,metabolite,replacement)
% replaceMets
%   Replaces metabolite names and annotation with replacement metabolite
%   that is already in the model. If this results in duplicate metabolites,
%   the replacement metabolite will be kept, while the S matrix is updated
%   to use the replacement metabolite instead.
%
%   model               a model structure
%   metabolite          string with name of metabolite to be replace
%   replacement         string with name of replacement metabolite
%
%   This function is useful when the model contains both 'oxygen' and 'o2'
%   as metabolites.
%
%   Usage: model=replaceMets(model,metabolite,replacement)
%
%   Eduard Kerkhoven, 2017-11-01

% Find occurence of replacement metabolites. Annotation will be taken from
% first metabolite found. Metabolite ID from replacement will be used where
% possible.
repIdx = find(strcmp(replacement,model.metNames));

% Change name and information from metabolite to replacement metabolite
metIdx = find(strcmp(metabolite,model.metNames));
model.metNames(metIdx) = model.metNames(repIdx(1));
if isfield(model,'metFormulas')
    model.metFormulas(metIdx) = model.metFormulas(repIdx(1));
end
if isfield(model,'metMiriams')
    model.metMiriams(metIdx) = model.metMiriams(repIdx(1));
end
if isfield(model,'metCharges')
    model.metCharges(metIdx) = model.metCharges(repIdx(1));
end
if isfield(model,'inchis')
    model.inchis(metIdx) = model.inchis(repIdx(1));
end
% Run through replacement metabolites and their compartments. If any of the
% to-be-replaced metabolites is already present (checked by
% metaboliteName[compartment], then the replacement metabolite is kept and
% the to-be-replace metabolite ID deleted.

% Build list of metaboliteName[compartment]
metCompsN =cellstr(num2str(model.metComps));
map = containers.Map(cellstr(num2str(transpose(1:length(model.comps)))),model.comps);
metCompsN = map.values(metCompsN);
metCompsN = strcat(lower(model.metNames),'[',metCompsN,']');

idxDelete=[];
for i = 1:length(repIdx)
    metCompsNidx=find(strcmp(metCompsN(repIdx(i)), metCompsN));
    if gt(length(metCompsNidx),1) % If more than 1 metabolite matches
        model.S(repIdx(i),:) = model.S(repIdx(i),:) + model.S(metCompsNidx(2:end),:);
        idxDelete=[idxDelete; metCompsNidx(2:end)]; % Make list of metabolite IDs to delete
    end
end

if ~isempty(idxDelete)
    model.S(idxDelete,:) =[];
    model.mets(idxDelete) = [];
    model.metNames(idxDelete) = [];
    model.metComps(idxDelete) = [];
    model.b(idxDelete) = [];
    if isfield(model,'metFormulas')
        model.metFormulas(idxDelete) = [];
    end
    if isfield(model,'unconstrained')
        model.unconstrained(idxDelete) = [];
    end
    if isfield(model,'metMiriams')
        model.metMiriams(idxDelete) = [];
    end
    if isfield(model,'metCharges')
        model.metCharges(idxDelete) = [];
    end
    if isfield(model,'inchis')
        model.inchis(idxDelete) = [];
    end
    if isfield(model,'metFrom')
        model.metFrom(idxDelete) = [];
    end
end

% This could now have created duplicate reactions. Contract model.
model=contractModel(model);
end