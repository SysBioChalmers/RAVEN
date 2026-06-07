function [solution, metabolite]=makeSomething(model,ignoreMets,isNames,minNrFluxes,allowExcretion,params,ignoreIntBounds)
% makeSomething
%   Tries to excrete any metabolite using as few reactions as possible.
%   The intended use is when you want to make sure that you model cannot
%   synthesize anything from nothing. It is then a faster way than to use
%   checkProduction or canProduce
%
%   model           a model structure
%   ignoreMets      either a cell array of metabolite IDs, a logical vector
%                   with the same number of elements as metabolites in the model,
%                   of a vector of indexes for metabolites to exclude from
%                   this analysis (optional, default [])
%   isNames         true if the supplied mets represent metabolite names
%                   (as opposed to IDs). This is a way to delete
%                   metabolites in several compartments at once without
%                   knowing the exact IDs. This only works if ignoreMets
%                   is a cell array (optional, default false)
%   minNrFluxes     solves the MILP problem of minimizing the number of
%                   fluxes instead of the sum. Slower, but can be
%                   used if the sum gives too many fluxes (optional, default
%                   false)
%   allowExcretion  allow for excretion of all other metabolites (optional,
%                   default true)
%   params          *obsolete option*
%   ignoreIntBounds	true if internal bounds (including reversibility)
%                   should be ignored. Exchange reactions are not affected.
%                   This can be used to find unbalanced solutions which are
%                   not possible using the default constraints (optional,
%                   default false)
%
%   solution        flux vector for the solution
%   metabolite      the index of the metabolite(s) which was excreted. If
%                   possible only one metabolite is reported, but there are
%                   situations where metabolites can only be excreted in
%                   pairs (or more)
%
%   NOTE: This works by forcing at least 1 unit of "any metabolites" to be
%   produced and then minimize for the sum of fluxes. If more than one
%   metabolite is produced, it picks one of them to be produced and then
%   minimizes for the sum of fluxes.
%
% Usage: [solution, metabolite]=makeSomething(model,ignoreMets,isNames,...
%           minNrFluxes,allowExcretion,params,ignoreIntBounds)

if nargin<2
    ignoreMets=[];
elseif ~islogical(ignoreMets) && ~isnumeric(ignoreMets)
    ignoreMets=convertCharArray(ignoreMets);
end
if nargin<3
    isNames=false;
end
if nargin<4
    minNrFluxes=false;
end
if nargin<5
    allowExcretion=true;
end
if nargin<6
    params.relGap=0.8;
end
if nargin<7
    ignoreIntBounds=false;
end

if isNames==true && ~isempty(ignoreMets)
    %Check that metsToRemove is a cell array
    if iscellstr(ignoreMets)==false
        EM='Must supply a cell array of strings if isNames=true';
        dispEM(EM);
    end
end

if isNames==false
    indexesToIgnore=getIndexes(model,ignoreMets,'mets');
else
    indexesToIgnore=[];
    for i=1:numel(ignoreMets)
        indexesToIgnore=[indexesToIgnore;find(strcmp(ignoreMets(i),model.metNames))];
    end
end

%Allow for excretion of all metabolites
if allowExcretion==true
    model.b=[model.b(:,1) inf(numel(model.mets),1)];
end

%Change all internal reactions to be unbounded in both directions
if ignoreIntBounds==true
    [~, I]=getExchangeRxns(model);
    nonExc=true(numel(model.rxns),1);
    nonExc(I)=false;
    model=setParam(model,'lb',nonExc,-1000);
    model=setParam(model,'ub',nonExc,1000);
    model=setParam(model,'rev',nonExc,1);
end

solution=[];
metabolite=[];

nRxns=numel(model.rxns);
nMets=numel(model.mets);

%Add excretion reactions for all metabolites
model.S=[model.S speye(nMets)*-1];

%Add so that they all produce a fake metabolite
model.S=[model.S;[sparse(1,nRxns) ones(1,nMets)]];

%Change so that the ignoreMets have a coefficient 0 in their reactions.
%Does not remove the actual reaction to make mapping easier later
model.S(:,indexesToIgnore+nRxns)=0;

%Add an excretion reaction for this last fake metabolite
model.S(size(model.S,1),size(model.S,2)+1)=-1;
model.b=[model.b;zeros(1,size(model.b,2))];
model.lb=[model.lb;zeros(nMets,1);1];
model.ub=[model.ub;inf(nMets+1,1)];
model.rev=[model.rev;zeros(nMets+1,1)];
model.c=zeros(size(model.S,2),1);

%Add padding to the reaction annotation to prevent an error in solveLP
padding=1:numel(model.rev);
padding=num2cell(padding)';
padding=cellfun(@num2str,padding,'uniformoutput',false);
model.rxns=padding;
model.rxnNames=padding;
model.eccodes=padding;
model.rxnMiriams=padding;
model.grRules=padding;
if isfield(model,'genes')
    model.rxnGeneMat=sparse(numel(model.rev),numel(model.genes));
end
model.subSystems=padding;
model.rxnFrom=padding;
model.rxnComps=ones(numel(model.rev),1);
model.rxnNotes=padding;
model.rxnReferences=padding;
model.rxnConfidenceScores=NaN(numel(model.rev),1);

sol=solveLP(model,1);
if any(sol.x)
    %It could be that several metabolites were produced in order to get the
    %best solution. The setdiff is to avoid including the last fake
    %metabolite
    I=setdiff(find(sol.x(nRxns+1:end)>0.1),size(model.S,1));
    
    if any(I) %This should always be true
        %Change the coefficients so that only the first is produced. This
        %is not always possible, but it is tested for since it it results
        %in more easily interpretable results
        
        oldS=model.S;
        foundSingle=false;
        %Test if any of the metabolites could be produced on their own
        for i=1:numel(I)
            model.S=oldS;
            J=nRxns+1:numel(model.lb)-1;
            J(I(i))=[];
            model.S(:,J)=0;
            sol=solveLP(model);
            if any(sol.x)
                foundSingle=true;
                break;
            end
        end
        %This means that there was no metabolite which could be produced on
        %its own. Then let all the produceable metabolites be produced.
        if foundSingle==false
            model.S=oldS;
        end
        if minNrFluxes==true
            %Has to add names for the rxns to prevent an error in
            %minNrFluxes
            model.rxns=cell(numel(model.lb),1);
            model.rxns(:)=sprintfc('tmp_%d',1:numel(model.lb));
            model.mets=cell(size(model.b,1),1);
            model.mets(:)={'TEMP'};
            sol=solveLP(model,3,params);
        else
            sol=solveLP(model,1);
        end
        solution=sol.x(1:nRxns);
        metabolite=find(sol.x(nRxns+1:end-1)>0.1);
    end
end
end
