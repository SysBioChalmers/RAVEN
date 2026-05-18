function model = initYAMLmodel(modelFields)
% initYAMLmodel
%   Build an empty RAVEN model struct from a Nx2 (fieldName, defaultValue)
%   schema (see legacyYAMLmodelFields). Used by both readYAMLmodel paths
%   before the parser fills the fields in.
model = struct();
for i=1:size(modelFields,1)
    model.(modelFields{i,1})=modelFields{i,2};
end
end
