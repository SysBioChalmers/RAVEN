function checkModelStruct(model,throwErrors,trimWarnings)
% checkModelStruct
%   Performs a number of checks to ensure that a model structure is ok
%
%   model           a model structure
%   throwErrors     true if the function should throw errors if
%                   inconsistencies are found. The alternative is to
%                   print warnings for all types of issues (opt, default true)
%   trimWarnings    true if only a maximal of 10 items should be displayed in
%                   a given error/warning (opt, default true)
%
%   NOTE: This is performed after importing a model from Excel or before
%   attempting to export a model to SBML format.
%
%   Usage: checkModelStruct(model,throwErrors,trimWarnings)
%
%   Simonas Marcisauskas, 2018-03-18
%

if nargin<2
    throwErrors=true;
end
if nargin<3
    trimWarnings=true;
end

%Missing elements
fields={'id';'description';'rxns';'mets';'S';'lb';'ub';'rev';'c';'b';'comps';'metComps'};
for i=1:numel(fields)
    if ~isfield(model,fields{i})
        EM=['The model is missing the "' fields{i} '" field'];
        dispEM(EM,throwErrors);
    end
end

%Type check
if ~ischar(model.id)
    EM='The "id" field must be a string';
    dispEM(EM,throwErrors);
end
if ~ischar(model.description)
    EM='The "description" field must be a string';
    dispEM(EM,throwErrors);
end
if ~iscellstr(model.rxns)
    EM='The "rxns" field must be a cell array of strings';
    dispEM(EM,throwErrors);
end
if ~iscellstr(model.mets)
    EM='The "mets" field must be a cell array of strings';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.S)
    EM='The "S" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.lb)
    EM='The "lb" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.ub)
    EM='The "ub" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.rev)
    EM='The "rev" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.c)
    EM='The "c" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.b)
    EM='The "b" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if ~iscellstr(model.comps)
    EM='The "comps" field must be a cell array of strings';
    dispEM(EM,throwErrors);
end
if ~isnumeric(model.metComps)
    EM='The "metComps" field must be of type "double"';
    dispEM(EM,throwErrors);
end
if isfield(model,'compNames')
    if ~iscellstr(model.compNames)
        EM='The "compNames" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'compOutside')
    if ~iscellstr(model.compOutside)
        EM='The "compOutside" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnNames')
    if ~iscellstr(model.rxnNames)
        EM='The "rxnNames" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'metNames')
    if ~iscellstr(model.metNames)
        EM='The "metNames" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'genes')
    if ~iscellstr(model.genes)
        EM='The "genes" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnGeneMat')
    if ~isnumeric(model.rxnGeneMat)
        EM='The "rxnGeneMat" field must be of type "double"';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'grRules')
    if ~iscellstr(model.grRules)
        EM='The "grRules" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnComps')
    if ~isnumeric(model.rxnComps)
        EM='The "rxnComps" field must be of type "double"';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'inchis')
    if ~iscellstr(model.inchis)
        EM='The "inchis" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'metFormulas')
    if ~iscellstr(model.metFormulas)
        EM='The "metFormulas" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'metCharges')
    if ~isnumeric(model.metCharges)
        EM='The "metCharges" field must be a double';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'subSystems')
    for i=1:numel(model.subSystems)
        if ~iscell(model.subSystems{i,1})
            EM='The "subSystems" field must be a cell array';
            dispEM(EM,throwErrors);
        end
    end
end
if isfield(model,'eccodes')
    if ~iscellstr(model.eccodes)
        EM='The "eccodes" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'unconstrained')
    if ~isnumeric(model.unconstrained)
        EM='The "unconstrained" field must be of type "double"';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnNotes')
    if ~iscellstr(model.rxnNotes)
        EM='The "rxnNotes" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnReferences')
    if ~iscellstr(model.rxnReferences)
        EM='The "rxnReferences" field must be a cell array of strings';
        dispEM(EM,throwErrors);
    end
end
if isfield(model,'rxnConfidenceScores')
    if ~isnumeric(model.rxnConfidenceScores)
        EM='The "rxnConfidenceScores" field must be a double';
        dispEM(EM,throwErrors);
    end
end

%Empty strings
if isempty(model.id)
    EM='The "id" field cannot be empty';
    dispEM(EM,throwErrors);
end
if any(cellfun(@isempty,model.rxns))
    EM='The model contains empty reaction IDs';
    dispEM(EM,throwErrors);
end
if any(cellfun(@isempty,model.mets))
    EM='The model contains empty metabolite IDs';
    dispEM(EM,throwErrors);
end
if any(cellfun(@isempty,model.comps))
    EM='The model contains empty compartment IDs';
    dispEM(EM,throwErrors);
end
EM='The following metabolites have empty names:';
dispEM(EM,throwErrors,model.mets(cellfun(@isempty,model.metNames)),trimWarnings);

if isfield(model,'genes')
    if any(cellfun(@isempty,model.genes))
        EM='The model contains empty gene IDs';
        dispEM(EM,throwErrors);
    end
end

%Duplicates
EM='The following reaction IDs are duplicates:';
dispEM(EM,throwErrors,model.rxns(duplicates(model.rxns)),trimWarnings);
EM='The following metabolite IDs are duplicates:';
dispEM(EM,throwErrors,model.mets(duplicates(model.mets)),trimWarnings);
EM='The following compartment IDs are duplicates:';
dispEM(EM,throwErrors,model.comps(duplicates(model.comps)),trimWarnings);
if isfield(model,'genes')
    EM='The following genes are duplicates:';
    dispEM(EM,throwErrors,model.genes(duplicates(model.genes)),trimWarnings);
end
metInComp=strcat(model.metNames,'[',model.comps(model.metComps),']');
EM='The following metabolites already exist in the same compartment:';
dispEM(EM,throwErrors,metInComp(duplicates(metInComp)),trimWarnings);

%Elements never used (print only as warnings
EM='The following reactions are empty (no involved metabolites):';
dispEM(EM,false,model.rxns(~any(model.S,1)),trimWarnings);
'The following metabolites are never used in a reaction:';
dispEM(EM,false,model.mets(~any(model.S,2)),trimWarnings);
if isfield(model,'genes')
    EM='The following genes are not associated to a reaction:';
    dispEM(EM,false,model.genes(~any(model.rxnGeneMat,1)),trimWarnings);
end
I=true(numel(model.comps),1);
I(model.metComps)=false;
EM='The following compartments contain no metabolites:';
dispEM(EM,false,model.comps(I),trimWarnings);

%Contradicting bounds
EM='The following reactions have contradicting bounds:';
dispEM(EM,throwErrors,model.rxns(model.lb>model.ub),trimWarnings);
EM='The following reactions have bounds contradicting their reversibility:';
dispEM(EM,throwErrors,model.rxns(model.lb<0 & model.rev==0),trimWarnings);

%Mapping of compartments
if isfield(model,'compOutside')
    EM='The following compartments are in "compOutside" but not in "comps":';
    dispEM(EM,throwErrors,setdiff(model.compOutside,[{''};model.comps]),trimWarnings);
end

%Met names which start with number
I=false(numel(model.metNames),1);
for i=1:numel(model.metNames)
    index=strfind(model.metNames{i},' ');
    if any(index)
        if any(str2double(model.metNames{i}(1:index(1)-1)))
            I(i)=true;
        end
    end
end
EM='The following metabolite names begin with a number directly followed by space:';
dispEM(EM,throwErrors,model.mets(I),trimWarnings);

%Non-parseable composition
if isfield(model,'metFormulas')
    [~, ~, exitFlag]=parseFormulas(model.metFormulas,true,false);
    EM='The composition for the following metabolites could not be parsed:';
    dispEM(EM,false,model.mets(exitFlag==-1),trimWarnings);
end

%Check if there are metabolites with different names but the same MIRIAM
%codes
if isfield(model,'metMiriams')
    miriams=containers.Map();
    for i=1:numel(model.mets)
        if ~isempty(model.metMiriams{i})
            %Loop through and add for each miriam
            for j=1:numel(model.metMiriams{i}.name)
                %Get existing metabolite indexes
                current=strcat(model.metMiriams{i}.name{j},'/',model.metMiriams{i}.value{j});
                if isKey(miriams,current)
                    existing=miriams(current);
                else
                    existing=[];
                end
                miriams(current)=[existing;i];
            end
        end
    end
    
    %Get all keys
    allMiriams=keys(miriams);
    
    hasMultiple=false(numel(allMiriams),1);
    for i=1:numel(allMiriams)
        if numel(miriams(allMiriams{i}))>1
            %Check if they all have the same name
            if numel(unique(model.metNames(miriams(allMiriams{i}))))>1
                hasMultiple(i)=true;
            end
        end
    end
    
    %Print output
    EM='The following MIRIAM strings are associated to more than one unique metabolite name:';
    dispEM(EM,false,allMiriams(hasMultiple));
end

%Check if there are metabolites with different names but the same InChI
%codes
if isfield(model,'inchis')
    inchis=containers.Map();
    for i=1:numel(model.mets)
        if ~isempty(model.inchis{i})
            %Get existing metabolite indexes
            if isKey(inchis,model.inchis{i})
                existing=inchis(model.inchis{i});
            else
                existing=[];
            end
            inchis(model.inchis{i})=[existing;i];
        end
    end
    
    %Get all keys
    allInchis=keys(inchis);
    
    hasMultiple=false(numel(allInchis),1);
    for i=1:numel(allInchis)
        if numel(inchis(allInchis{i}))>1
            %Check if they all have the same name
            if numel(unique(model.metNames(inchis(allInchis{i}))))>1
                hasMultiple(i)=true;
            end
        end
    end
    
    %Print output
    EM='The following InChI strings are associated to more than one unique metabolite name:';
    dispEM(EM,false,allInchis(hasMultiple));
end
end

function I=duplicates(strings)
I=false(numel(strings),1);
[J, K]=unique(strings);
if numel(J)~=numel(strings)
    L=1:numel(strings);
    L(K)=[];
    I(L)=true;
end
end
