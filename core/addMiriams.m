function model=addMiriams(model,type,id,name,value)
% addMiriams
%   Adds MIRIAM annotations to an existing model.
%
%   model           model MIRIAMs should be added
%   type            'rxns' or 'mets' for reaction or metabolite annotations
%   id              character string or string array with reaction or
%                   metabolite ids to which MIRIAMs should be written
%   name            character string with name of the source of MIRIAM
%                   annotation, where possible compliant with
%                   identifiers.org (e.g. kegg.compound)
%   value           character string or string array with identifiers.
%                   If 'id' is a string array, then 'value' should be of
%                   the same length (assuming same order), or 'value'
%                   should be a character string (assuming same 'value' for
%                   all specified ids).
%
%   model           model structure with updated MIRIAM annotations
%
%   Annotations can be provided for one identifier-name at the time (e.g.
%   kegg.compound, chebi or sbo). However, values can be provided in
%   various ways:
%   -   one annotation for one metabolite or reaction: id and value are
%       both character strings
%   -   identical annotation for multiple mets or rxns: id is string array
%       while value is character string
%   -   different annotations for multiple mets or rxns: id and value are
%       both string arrays
%
%   Usage: model=addMiriams(model,type,id,name,value)
%
%   Eduard Kerkhoven, 2018-09-15

% Check if valid input is provided
if ~ismember(type,{'mets','rxns'})
    error('type should be specified as ''mets'' or ''rxns''')
end

% If id is array, but value is string, then assume same value for all id
if iscellstr(id) && ischar(value)
    value2(1:length(id),1)={value};
    value=value2; clear value2
end

% If id is char, convert to array
if ischar(id)
    id={id};
end

% Convert id to indexes
if type=='mets'
    id = getIndexes(model,id,'mets');
else
    id = getIndexes(model,id,'rxns');
end
for i=1:length(id)
    if type=='mets'
        Miriam=model.metMiriams{i};
    else
        Miriam=model.rxnMiriams{i};
    end
    miriamIdx=ismember(Miriam.name,name);
    % Add annotation only if it doesn't exist yet
    if ~any(miriamIdx) || ~ismember(Miriam.value(miriamIdx),value(i))
        Miriam.name(end+1)={name};
        Miriam.value(end+1)=value(i);
    end
    if type=='mets'
        model.metMiriams{i}=Miriam;
    else
        model.rxnMiriams{i}=Miriam;
    end
end
end
