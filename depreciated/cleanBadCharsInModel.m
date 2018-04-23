function model=cleanBadCharsInModel(model)
% cleanBadCharsInModel
%   Converts characters which are illegal in SBML to their entity reference
%
%   model	a model structure
%
%   model   a model structure without characters which are illegal in SBML
%
%   The characters <>&'" are illegal in SBML (in most cases). They are here
%   converted to their entity references. The main purpose of this function
%   is for using in exportModel to ensure that the resulting SBML file is
%   valid.
%
%   Usage: model=cleanBadCharsInModel(model)
%
%   Rasmus Agren, 2013-08-03
%

%These are IDs and shouldn't have to be checked since it's illegal to have
%non-alphanumeric characters there. This should be checked in
%checkModelStruct instead.
model.id=cleanBadChars(model.id);
model.mets=cleanBadChars(model.mets);
model.rxns=cleanBadChars(model.rxns);
model.comps=cleanBadChars(model.comps);
if isfield(model,'compOutside')
    model.compOutside=cleanBadChars(model.compOutside);
end

if isfield(model,'compNames')
    model.compNames=cleanBadChars(model.compNames);
end
if isfield(model,'description')
    model.description=cleanBadChars(model.description);
end
if isfield(model,'rxnNames')
    model.rxnNames=cleanBadChars(model.rxnNames);
end
if isfield(model,'grRules')
    model.grRules=cleanBadChars(model.grRules);
end
if isfield(model,'subSystems')
    model.subSystems=cleanBadChars(model.subSystems);
end
if isfield(model,'eccodes')
    model.eccodes=cleanBadChars(model.eccodes);
end
if isfield(model,'metNames')
    model.metNames=cleanBadChars(model.metNames);
end
if isfield(model,'genes')
    model.genes=cleanBadChars(model.genes);
end
if isfield(model,'geneShortName')
    model.geneShortName=cleanBadChars(model.geneShortName);
end
if isfield(model,'inchis')
    model.inchis=cleanBadChars(model.inchis);
end
if isfield(model,'metFormulas')
    model.metFormulas=cleanBadChars(model.metFormulas);
end

%These fields are cell arrays of structures
if isfield(model,'metMiriams')
    for i=1:numel(model.metMiriams)
        if ~isempty(model.metMiriams{i})
            model.metMiriams{i}.name=cleanBadChars(model.metMiriams{i}.name);
            model.metMiriams{i}.value=cleanBadChars(model.metMiriams{i}.value);
        end
    end
end
if isfield(model,'rxnMiriams')
    for i=1:numel(model.rxnMiriams)
        if ~isempty(model.rxnMiriams{i})
            model.rxnMiriams{i}.name=cleanBadChars(model.rxnMiriams{i}.name);
            model.rxnMiriams{i}.value=cleanBadChars(model.rxnMiriams{i}.value);
        end
    end
end
if isfield(model,'compMiriams')
    for i=1:numel(model.compMiriams)
        if ~isempty(model.compMiriams{i})
            model.compMiriams{i}.name=cleanBadChars(model.compMiriams{i}.name);
            model.compMiriams{i}.value=cleanBadChars(model.compMiriams{i}.value);
        end
    end
end
if isfield(model,'geneMiriams')
    for i=1:numel(model.geneMiriams)
        if ~isempty(model.geneMiriams{i})
            model.geneMiriams{i}.name=cleanBadChars(model.geneMiriams{i}.name);
            model.geneMiriams{i}.value=cleanBadChars(model.geneMiriams{i}.value);
        end
    end
end
end
