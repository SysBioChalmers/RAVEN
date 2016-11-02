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
%   Rasmus Agren, 2013-11-03
%   Simonas Marcisauskas, 2016-11-01 - added checks for rxnNotes,
%   rxnReferences, confidenceScores and metCharge
%

if nargin<2
    throwErrors=true;
end
if nargin<3
    trimWarnings=true;
end

%Missing elements
if ~isfield(model,'id')
    dispEM('The model is missing the "id" field',throwErrors);
end
if ~isfield(model,'description')
    dispEM('The model is missing the "description" field',throwErrors);
end
if ~isfield(model,'rxns')
    dispEM('The model is missing the "rxns" field',throwErrors);
end
if ~isfield(model,'mets')
    dispEM('The model is missing the "mets" field',throwErrors);
end
if ~isfield(model,'S')
    dispEM('The model is missing the "S" field',throwErrors);
end
if ~isfield(model,'lb')
    dispEM('The model is missing the "lb" field',throwErrors);
end
if ~isfield(model,'ub')
    dispEM('The model is missing the "ub" field',throwErrors);
end
if ~isfield(model,'rev')
    dispEM('The model is missing the "rev" field',throwErrors);
end
if ~isfield(model,'c')
    dispEM('The model is missing the "c" field',throwErrors);
end
if ~isfield(model,'b')
    dispEM('The model is missing the "b" field',throwErrors);
end
if ~isfield(model,'comps')
    dispEM('The model is missing the "comps" field',throwErrors);
end
if ~isfield(model,'metComps')
    dispEM('The model is missing the "metComps" field',throwErrors);
end

%Type check
if ~ischar(model.id)
    dispEM('The "id" field must be a string',throwErrors);
end
if ~ischar(model.description)
    dispEM('The "description" field must be a string',throwErrors);
end
if ~iscellstr(model.rxns)
    dispEM('The "rxns" field must be a cell array of strings',throwErrors);
end
if ~iscellstr(model.mets)
    dispEM('The "mets" field must be a cell array of strings',throwErrors);
end
if ~isnumeric(model.S)
    dispEM('The "S" field must be of type "double"',throwErrors);
end
if ~isnumeric(model.lb)
    dispEM('The "lb" field must be of type "double"',throwErrors);
end
if ~isnumeric(model.ub)
    dispEM('The "ub" field must be of type "double"',throwErrors);
end
if ~isnumeric(model.rev)
    dispEM('The "rev" field must be of type "double"',throwErrors);
end
if ~isnumeric(model.c)
    dispEM('The "c" field must be of type "double"',throwErrors);
end
if ~isnumeric(model.b)
    dispEM('The "b" field must be of type "double"',throwErrors);
end
if ~iscellstr(model.comps)
    dispEM('The "comps" field must be a cell array of strings',throwErrors);
end
if ~isnumeric(model.metComps)
    dispEM('The "metComps" field must be of type "double"',throwErrors);
end
if isfield(model,'compNames')
    if ~iscellstr(model.compNames)
        dispEM('The "compNames" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'compOutside')
    if ~iscellstr(model.compOutside)
        dispEM('The "compOutside" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'rxnNames')
    if ~iscellstr(model.rxnNames)
        dispEM('The "rxnNames" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'metNames')
    if ~iscellstr(model.metNames)
        dispEM('The "metNames" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'genes')
    if ~iscellstr(model.genes)
        dispEM('The "genes" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'rxnGeneMat')
    if ~isnumeric(model.rxnGeneMat)
        dispEM('The "rxnGeneMat" field must be of type "double"',throwErrors);
    end
end
if isfield(model,'grRules')
    if ~iscellstr(model.grRules)
        dispEM('The "grRules" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'rxnComps')
    if ~isnumeric(model.rxnComps)
        dispEM('The "rxnComps" field must be of type "double"',throwErrors);
    end
end
if isfield(model,'inchis')
    if ~iscellstr(model.inchis)
        dispEM('The "inchis" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'metFormulas')
    if ~iscellstr(model.metFormulas)
        dispEM('The "metFormulas" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'metCharge')
    if ~isnumeric(model.metCharge)
        dispEM('The "metCharge" field must be of type "double"',throwErrors);
    end
end
if isfield(model,'subSystems')
    if ~iscellstr(model.subSystems)
        dispEM('The "subSystems" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'eccodes')
    if ~iscellstr(model.eccodes)
        dispEM('The "eccodes" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'unconstrained')
    if ~isnumeric(model.unconstrained)
        dispEM('The "unconstrained" field must be of type "double"',throwErrors);
    end
end
if isfield(model,'rxnNotes')
    if ~iscellstr(model.rxnNotes)
        dispEM('The "rxnNotes" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'rxnReferences')
    if ~iscellstr(model.rxnReferences)
        dispEM('The "rxnReferences" field must be a cell array of strings',throwErrors);
    end
end
if isfield(model,'confidenceScores')
    if ~iscellstr(model.confidenceScores)
        dispEM('The "confidenceScores" field must be a cell array of strings',throwErrors);
    end
end

%Empty strings
if isempty(model.id)
    dispEM('The "id" field cannot be empty',throwErrors);
end
if any(cellfun(@isempty,model.rxns))
    dispEM('The model contains empty reaction IDs',throwErrors);
end
if any(cellfun(@isempty,model.mets))
    dispEM('The model contains empty metabolite IDs',throwErrors);
end
if any(cellfun(@isempty,model.comps))
    dispEM('The model contains empty compartment IDs',throwErrors);
end
dispEM('The following metabolites have empty names:',throwErrors,model.mets(cellfun(@isempty,model.metNames)),trimWarnings);

if isfield(model,'genes')
    if any(cellfun(@isempty,model.genes))
        dispEM('The model contains empty gene IDs',throwErrors);
    end
end

%Illegal characters in IDs
dispEM('Illegal characters in reaction IDs:',throwErrors,model.rxns(illegal(model.rxns,'id')),trimWarnings);
dispEM('Illegal characters in metabolite IDs:',throwErrors,model.mets(illegal(model.mets,'id')),trimWarnings);
dispEM('Illegal characters in compartment IDs:',throwErrors,model.comps(illegal(model.comps,'id')),trimWarnings);

%Duplicates
dispEM('The following reaction IDs are duplicates:',throwErrors,model.rxns(duplicates(model.rxns)),trimWarnings);
dispEM('The following metabolite IDs are duplicates:',throwErrors,model.mets(duplicates(model.mets)),trimWarnings);
dispEM('The following compartment IDs are duplicates:',throwErrors,model.comps(duplicates(model.comps)),trimWarnings);
if isfield(model,'genes')
    dispEM('The following genes are duplicates:',throwErrors,model.genes(duplicates(model.genes)),trimWarnings);
end
metInComp=strcat(model.metNames,'[',model.comps(model.metComps),']');
dispEM('The following metabolites already exist in the same compartment:',throwErrors,metInComp(duplicates(metInComp)),trimWarnings);
    
%Elements never used (print only as warnings
dispEM('The following reactions are empty (no involved metabolites):',false,model.rxns(~any(model.S,1)),trimWarnings);
dispEM('The following metabolites are never used in a reaction:',false,model.mets(~any(model.S,2)),trimWarnings);
if isfield(model,'genes')
    dispEM('The following genes are not associated to a reaction:',false,model.genes(~any(model.rxnGeneMat,1)),trimWarnings);
end
I=true(numel(model.comps),1);
I(model.metComps)=false;
dispEM('The following compartments contain no metabolites:',false,model.comps(I),trimWarnings);

%Contradicting bounds
dispEM('The following reactions have contradicting bounds:',throwErrors,model.rxns(model.lb>model.ub),trimWarnings);
dispEM('The following reactions have bounds contradicting their reversibility:',throwErrors,model.rxns(model.lb<0 & model.rev==0),trimWarnings);

%Mapping of compartments
if isfield(model,'compOutside')
    dispEM('The following compartments are in "compOutside" but not in "comps":',throwErrors,setdiff(model.compOutside,[{''};model.comps]),trimWarnings);
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
dispEM('The following metabolite names begin with a number directly followed by space:',throwErrors,model.mets(I),trimWarnings);

%Non-parseable composition
if isfield(model,'metFormulas')
    [crap, crap, exitFlag]=parseFormulas(model.metFormulas,true,false);
     dispEM('The composition for the following metabolites could not be parsed:',false,model.mets(exitFlag==-1),trimWarnings);
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
             current=strcat(model.metMiriams{i}.name{j},':',model.metMiriams{i}.value{j});
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
    dispEM('The following MIRIAM strings are associated to more than one unique metabolite name:',false,allMiriams(hasMultiple));
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
    dispEM('The following InChI strings are associated to more than one unique metabolite name:',false,allInchis(hasMultiple));
end
end

function I=duplicates(strings)
    I=false(numel(strings),1);
    [J K]=unique(strings);
    if numel(J)~=numel(strings)
        L=1:numel(strings);
        L(K)=[];
        I(L)=true;
    end
end
function I=illegal(strings,type)
    %Just to save some space
    if strcmpi(type,'id')
        %Checks which strings in a cell array contains illegal characters
        I=cellfun(@any,regexp(strings,'[^a-z_A-Z0-9-]', 'once'));
    else

    end
end