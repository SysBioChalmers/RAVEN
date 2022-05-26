function newModel=addMets(model,metsToAdd,copyInfo,prefix)
% addMets
%   Adds metabolites to a model
%
%   model       a model structure
%   metsToAdd   the metabolite structure can have the following fields:
%               mets           cell array with unique strings that
%                              identifies each metabolite (opt, default IDs
%                              of new metabolites are numbered with the
%                              prefix defined below)
%               metNames       cell array with the names of each
%                              metabolite
%               compartments   cell array with the compartment of each
%                              metabolite. Should match model.comps.
%                              If this is a string rather than a cell
%                              array it is assumed that all mets are in
%                              that compartment
%               b              Nx1 or Nx2 matrix with equality constraints
%                              for each metabolite (opt, default 0)
%               unconstrained  vector describing if each metabolite is an
%                              exchange metabolite (1) or not (0) (opt,
%                              default 0)
%               inchis         cell array with InChI strings for each
%                              metabolite (opt, default '')
%               metFormulas    cell array with the formulas for each of
%                              the metabolites (opt, default '')
%               metMiriams     cell array with MIRIAM structures (opt,
%                              default [])
%               metCharges     metabolite charge (opt, default NaN)
%               metNotes       cell array with metabolite notes as strings
%                              (opt, default '')
%   copyInfo    when adding metabolites to a compartment where it previously
%               doesn't exist, the function will copy any available annotation
%               from the metabolite in another compartment (opt, default true)
%   prefix      when metsToAdd.mets is not specified, new metabolite IDs
%               are generated with the prefix specified here. If IDs with
%               the prefix are already used in the model then the
%               numbering will start from the highest existing integer+1
%               (opt, default 'm_')
%
%   newModel     an updated model structure
%
%   NOTE: This function does not make extensive checks about MIRIAM formats,
%   forbidden characters or such.
%
%   Example: If multiple metabolites are added at once, the metMiriams cell
%   array should be defined as (example with ChEBI and KEGG):
%
%   metsToAdd.metMiriams{1} = struct('name',{{'chebi';'kegg.compound'}},...
%       'value',{{'CHEBI:18072';'C11821'}});
%   metsToAdd.metMiriams{2} = struct('name',{{'chebi';'kegg.compound'}},...
%       'value',{{'CHEBI:31132';'C12248'}});
%
%   Usage: newModel=addMets(model,metsToAdd,copyInfo,prefix)

if nargin<3
    copyInfo=true;
end
if nargin<4
    prefix='m_';
else
    prefix=char(prefix);
end

newModel=model;

if isempty(metsToAdd)
    return;
end

%Check some stuff regarding the required fields
if ~isfield(metsToAdd,'mets')
    metsToAdd.mets=generateNewIds(newModel,'mets',prefix,numel(metsToAdd.metNames));
else
    metsToAdd.mets=convertCharArray(metsToAdd.mets);
end
if ~isfield(metsToAdd,'metNames')
    metsToAdd.metNames=metsToAdd.mets;
else
    metsToAdd.metNames=convertCharArray(metsToAdd.metNames);
end
if ~isfield(metsToAdd,'compartments')
    EM='compartments is a required field in metsToAdd';
    dispEM(EM);
else
    metsToAdd.compartments=convertCharArray(metsToAdd.compartments);
    %If only one compartment is given, assume it is for all metabolites
    if numel(metsToAdd.compartments)==1 && numel(metsToAdd.mets)>1
        temp=cell(numel(metsToAdd.mets),1);
        temp(:)={metsToAdd.compartments};
        metsToAdd.compartments=temp;
    end
end

%Number of metabolites
nMets=numel(metsToAdd.mets);
nOldMets=numel(model.mets);
filler=cell(nMets,1);
filler(:)={''};
largeFiller=cell(nOldMets,1);
largeFiller(:)={''};

%Check that no metabolite ids are already present in the model
I=ismember(metsToAdd.mets,model.mets);
if any(I)
    dispEM(['One or more elements in metsToAdd.mets are already present in model.mets: '...
        metsToAdd.mets{I}]);
else
    newModel.mets=[newModel.mets;metsToAdd.mets(:)];
end

%Check that all the compartments could be found
[I, compMap]=ismember(metsToAdd.compartments,model.comps);
if ~all(I)
    EM='metsToAdd.compartments must match model.comps';
    dispEM(EM);
end

%Check that the metabolite names aren't present in the same compartment.
%Not the neatest way maybe..
t1=strcat(metsToAdd.metNames(:),'***',metsToAdd.compartments(:));
t2=strcat(model.metNames,'***',model.comps(model.metComps));
if any(ismember(t1,t2))
    EM='One or more elements in metsToAdd.metNames already exist in the same compartments as the one it is being added to';
    dispEM(EM);
end

%Some more checks and if they pass then add each field to the structure
if numel(metsToAdd.metNames)~=nMets
    EM='metsToAdd.metNames must have the same number of elements as metsToAdd.mets';
    dispEM(EM);
else
    newModel.metNames=[newModel.metNames;metsToAdd.metNames(:)];
end

if numel(compMap)~=nMets
    EM='metsToAdd.compartments must have the same number of elements as metsToAdd.mets';
    dispEM(EM);
else
    newModel.metComps=[newModel.metComps;compMap(:)];
end

if isfield(metsToAdd,'b')
    if size(metsToAdd.b,1)~=nMets
        EM='metsToAdd.b must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    else
        %Add empty field if it doesn't exist
        if ~isfield(newModel,'b')
            newModel.b=zeros(nOldMets,1);
        end
        
        %If the original is only one vector
        if size(metsToAdd.b,2)>size(newModel.b,2)
            newModel.b=[newModel.b newModel.b];
        end
        %Add the new ones
        newModel.b=[newModel.b;metsToAdd.b];
    end
else
    if isfield(newModel,'b')
        %Add the default
        newModel.b=[newModel.b;zeros(nMets,size(newModel.b,2))];
    end
end

if isfield(metsToAdd,'unconstrained')
    if numel(metsToAdd.unconstrained)~=nMets
        EM='metsToAdd.unconstrained must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    else
        %Add empty field if it doesn't exist
        if ~isfield(newModel,'unconstrained')
            newModel.unconstrained=zeros(nOldMets,1);
        end
        
        %Add the new ones
        newModel.unconstrained=[newModel.unconstrained;metsToAdd.unconstrained(:)];
    end
else
    if isfield(newModel,'unconstrained')
        %Add the default
        newModel.unconstrained=[newModel.unconstrained;zeros(nMets,1)];
    end
end

if isfield(metsToAdd,'inchis')
    metsToAdd.inchis=convertCharArray(metsToAdd.inchis);
    if numel(metsToAdd.inchis)~=nMets
        EM='metsToAdd.inchis must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'inchis')
        newModel.inchis=largeFiller;
    end
    newModel.inchis=[newModel.inchis;metsToAdd.inchis(:)];
else
    %Add empty strings if structure is in model
    if isfield(newModel,'inchis')
        newModel.inchis=[newModel.inchis;filler];
    end
end

if isfield(metsToAdd,'metFormulas')
    metsToAdd.metFormulas=convertCharArray(metsToAdd.metFormulas);
    if numel(metsToAdd.metFormulas)~=nMets
        EM='metsToAdd.metFormulas must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'metFormulas')
        newModel.metFormulas=largeFiller;
    end
    newModel.metFormulas=[newModel.metFormulas;metsToAdd.metFormulas(:)];
else
    %Add default
    if isfield(newModel,'metFormulas')
        newModel.metFormulas=[newModel.metFormulas;filler];
    end
end

if isfield(metsToAdd,'metCharges')
    if numel(metsToAdd.metCharges)~=nMets
        EM='metsToAdd.metCharges must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    end
    if ~isnumeric(metsToAdd.metCharges)
        EM='metsToAdd.metCharges must be of type "double"';
        dispEM(EM);
    end
    if ~isfield(newModel,'metCharges')
        newModel.metCharges=NaN(numel(largeFiller),1);
    end
    newModel.metCharges=[newModel.metCharges;metsToAdd.metCharges(:)];
else
    %Add default
    if isfield(newModel,'metCharges')
        newModel.metCharges=[newModel.metCharges;NaN(numel(filler),1)];
    end
end

if isfield(metsToAdd,'metNotes')
    metsToAdd.metNotes=convertCharArray(metsToAdd.metNotes);
    if numel(metsToAdd.metNotes)~=nMets
        EM='metsToAdd.metNotes must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'metNotes')
        newModel.metNotes=largeFiller;
    end
    newModel.metNotes=[newModel.metNotes;metsToAdd.metNotes(:)];
else
    %Add empty strings if structure is in model
    if isfield(newModel,'metNotes')
        newModel.metNotes=[newModel.metNotes;filler];
    end
end

%Don't check the type of metMiriams
if isfield(metsToAdd,'metMiriams')
    if numel(metsToAdd.metMiriams)~=nMets
        EM='metsToAdd.metMiriams must have the same number of elements as metsToAdd.mets';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'metMiriams')
        newModel.metMiriams=cell(nOldMets,1);
    end
    newModel.metMiriams=[newModel.metMiriams;metsToAdd.metMiriams(:)];
else
    if isfield(newModel,'metMiriams')
        newModel.metMiriams=[newModel.metMiriams;cell(nMets,1)];
    end
end

if isfield(newModel,'metFrom')
    newModel.metFrom=[newModel.metFrom;filler];
end

%Expand the S matrix
newModel.S=[newModel.S;sparse(nMets,size(newModel.S,2))];

if copyInfo==true
    [I, J]=ismember(metsToAdd.metNames,model.metNames);
    J=J(I);
    %I is the indexes of the new metabolites for which a metabolite with
    %the same name existed
    I=find(I)+nOldMets;
    %Go through each of the added mets and copy annotation if it doesn't
    %exist
    for i=1:numel(I)
        if isfield(newModel,'inchis')
            if isempty(newModel.inchis{I(i)})
                newModel.inchis(I(i))=newModel.inchis(J(i));
            end
        end
        if isfield(newModel,'metFormulas')
            if isempty(newModel.metFormulas{I(i)})
                newModel.metFormulas(I(i))=newModel.metFormulas(J(i));
            end
        end
        if isfield(newModel,'metMiriams')
            if isempty(newModel.metMiriams{I(i)})
                newModel.metMiriams(I(i))=newModel.metMiriams(J(i));
            end
        end
        if isfield(newModel,'metCharges')
            newModel.metCharges(I(i))=newModel.metCharges(J(i));
        end
    end
end
end
