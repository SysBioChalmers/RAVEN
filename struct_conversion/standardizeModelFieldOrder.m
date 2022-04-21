function orderedModel=standardizeModelFieldOrder(model)
% standardizeModelFieldOrder
%   Orders fields of RAVEN model structure as specified at
%   https://github.com/SysBioChalmers/RAVEN/wiki/RAVEN-Model-Structure
%
%   Input: model           model structure, either RAVEN or COBRA format
%
%   Output: orderedModel   model structure with ordered fields
%
%   The model fields themselves are not changed, only the order is
%   modified. For changing model fields between RAVEN and COBRA format, use
%   ravenCobraWrapper().
%
%   Usage: orderedModel=standardizeModelFieldOrder(model)

ravenPath=findRAVENroot();

if ~isfield(model,'rules') % Check if model is RAVEN
    fid = fopen(fullfile(ravenPath,'struct_conversion','orderRavenFields.csv'));
    fields = textscan(fid,'%s','Delimiter',',','HeaderLines',0);
    fields = fields{1};
    fclose(fid);
else % If model is COBRA
    fid = fopen(fullfile(ravenPath,'struct_conversion','COBRA_structure_fields.csv')); % Taken from https://github.com/opencobra/cobratoolbox/blob/develop/src/base/io/definitions/COBRA_structure_fields.csv
    fields = textscan(fid,repmat('%s',1,15),'Delimiter','\t','HeaderLines',1);
    fields = fields{1};
    fclose(fid);
end

modelfields = fieldnames(model);
order = fields(ismember(fields(:,1),modelfields));
remainingOrder = sort(setdiff(modelfields,order));
overallOrder = [columnVector(order);columnVector(remainingOrder)];
orderedModel = orderfields(model,overallOrder);
end