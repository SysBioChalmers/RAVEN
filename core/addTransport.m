function [model, addedRxns]=addTransport(model,fromComp,toComps,metNames,isRev,onlyToExisting,prefix)
% addTransport
%   Adds transport reactions between compartments
%
%   model           a model structure
%   fromComp        the id of the compartment to transport from (should
%                   match model.comps)
%   toComps         a cell array of compartment names to transport to (should
%                   match model.comps)
%   metNames        the metabolite names to add transport for (opt, all
%                   metabolites in fromComp)
%   isRev           true if the transport reactions should be reversible
%                   (opt, default true)
%   onlyToExisting  true if transport of a metabolite should only be added
%                   if it already exists in toComp. If false, then new metabolites
%                   are added with addMets first (opt, default true)
%   prefix          string specifying prefix to reaction IDs (opt, default
%                   'tr_')
%
%   This is a faster version than addRxns when adding transport reactions.
%   New reaction names are formatted as "metaboliteName, fromComp-toComp", 
%   while new reaction IDs are sequentially counted with a tr_ prefix:
%   e.g. tr_0001, tr_0002, etc.
%
%   Usage: [model, addedRxns]=addTransport(model,fromComp,toComps,metNames,...
%           isRev,onlyToExisting,prefix)

fromComp=char(fromComp);
[I, fromID]=ismember(model.comps,fromComp);
fromID=find(fromID);
if sum(I)~=1
    EM='fromComps must have exactly one match in model.comps';
    dispEM(EM);
end
toComps=convertCharArray(toComps);
[I, toIDs]=ismember(toComps,model.comps);
if ~all(I)
    EM='All compartments in toComps must have a match in model.comps';
    dispEM(EM);
end
if nargin<4
    %Find all metabolites in fromComp
    metNames=model.metNames(model.metComps==fromID);
elseif isempty(metNames)
    %Find all metabolites in fromComp
    metNames=model.metNames(ismember(model.metComps,model.comps(fromID)));
else
    metNames=convertCharArray(metNames);
end

if nargin<5
    isRev=true;
end
if nargin<6
    onlyToExisting=true;
end
if nargin<7
    prefix='tr_';
else
    prefix=char(prefix);
end

%Check that the names are unique
if numel(unique(metNames))~=numel(metNames)
    dispEM('Not all metabolite names are unique');
end

%Get the indexes of the mets in fromComp
I=find(model.metComps==fromID);
[J, K]=ismember(metNames,model.metNames(I));
if ~all(J)
    EM='Not all metabolites in metNames exist in fromComp';
    dispEM(EM);
end
fromMets=I(K); %These are the ids of the metabolites to transport. The order corresponds to metNames
addedRxns={};
%Loop through and add for each compartment in toComps
for i=1:numel(toComps)
    fromMetsInComp=fromMets; %If onlyToExisting==true then not all mets are transported to each compartment
    %Get the indexes of the mets in the compartment
    I=find(model.metComps==toIDs(i));
    [J, K]=ismember(metNames,model.metNames(I));
    if onlyToExisting==true || all(J)
        toMets=I(K(J)); %Only look at the existing ones
        fromMetsInComp=fromMetsInComp(J);
    else
        %This is if not all metabolites exist in the target compartment,
        %and they should be added
        metsToAdd.metNames=metNames(J==0);
        metsToAdd.compartments=toComps{i};
        model=addMets(model,metsToAdd);
        
        %Redo the mapping when all mets are there. A bit lazy, but it's
        %fast anyways
        I=find(model.metComps==toIDs(i));
        [~, K]=ismember(metNames,model.metNames(I));
        toMets=I(K); %All are guaranteed to be found now
    end
    
    %Construct the S matrix
    nRxns=numel(fromMetsInComp);
    newS=zeros(numel(model.mets),nRxns);
    newS(sub2ind(size(newS),fromMetsInComp(:),(1:nRxns)'))=-1;
    newS(sub2ind(size(newS),toMets(:),(1:nRxns)'))=1;
    
    %Add the reactions
    model.S=[model.S sparse(newS)];
    if isfield(model.annotation,'defaultLB')
        lb = model.annotation.defaultLB;
        ub = model.annotation.defaultUB;
    else
        lb = -inf;
        ub = inf;
    end
    if isRev==true
        model.lb=[model.lb;ones(nRxns,1)*lb];
        model.rev=[model.rev;ones(nRxns,1)];
    else
        model.lb=[model.lb;zeros(nRxns,1)];
        model.rev=[model.rev;zeros(nRxns,1)];
    end
    model.ub=[model.ub;ones(nRxns,1)*ub];
    model.c=[model.c;zeros(nRxns,1)];
    
    %Add annotation
    filler=cell(nRxns,1);
    filler(:)={''};
    addedRxnsID=generateNewIds(model,'rxns',prefix,nRxns);
    addedRxnsName=transpose(strcat(metNames, {' transport, '}, model.compNames(fromID), '-', model.compNames(toIDs(i))));
    model.rxns=[model.rxns;addedRxnsID];
    model.rxnNames=[model.rxnNames;addedRxnsName];
    
    if isfield(model,'eccodes')
        model.eccodes=[model.eccodes;filler];
    end
    if isfield(model,'subSystems')
        ssFiller=filler;
        if isRev==1
            ssFiller(:)={{['Transport between ' fromComp ' and ' toComps{i}]}};
        else
            ssFiller(:)={{['Transport from ' fromComp ' to ' toComps{i}]}};
        end
        model.subSystems=[model.subSystems;ssFiller];
    end
    if isfield(model,'grRules')
        model.grRules=[model.grRules;filler];
    end
    if isfield(model,'rxnFrom')
        model.rxnFrom=[model.rxnFrom;filler];
    end
    if isfield(model,'rxnMiriams')
        model.rxnMiriams=[model.rxnMiriams;cell(nRxns,1)];
    end
    if isfield(model,'rxnComps')
        model.rxnComps=[model.rxnComps;ones(nRxns,1)];
        fprintf('NOTE: The added transport reactions will be added to the first compartment\n');
    end
    if isfield(model,'rxnGeneMat')
        model.rxnGeneMat=[model.rxnGeneMat;sparse(nRxns,numel(model.genes))];
    end
    if isfield(model,'rxnNotes')
        model.rxnNotes=[model.rxnNotes;filler];
    end
    if isfield(model,'rxnReferences')
        model.rxnReferences=[model.rxnReferences;filler];
    end
    if isfield(model,'rxnConfidenceScores')
        model.rxnConfidenceScores=[model.rxnConfidenceScores;ones(nRxns,1)];
    end
    if isfield(model,'rxnDeltaG')
        model.rxnDeltaG=[model.rxnDeltaG;zeros(nRxns,1)];
    end
    addedRxns = [addedRxns; addedRxnsID];
end
end
