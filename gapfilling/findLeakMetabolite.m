function [solution, metabolite]=findLeakMetabolite(model, direction, varargin)
% findLeakMetabolite  Find a metabolite that can be freely made or consumed.
%
% Tests whether the model can produce (excrete) or consume (take up) any
% metabolite without balanced stoichiometric constraints — a sign that the
% model contains a stoichiometric leak, futile cycle, or unconstrained
% exchange. This is the unified replacement for makeSomething ('produce')
% and consumeSomething ('consume').
%
% Parameters
% ----------
% model : struct
%     a model structure.
% direction : char
%     'produce' to look for freely excreted metabolites, or 'consume' to
%     look for freely consumed metabolites.
%
% Name-Value Arguments
% --------------------
% ignoreMets : cell or logical or double
%     either a cell array of metabolite IDs, a logical vector with the same
%     number of elements as metabolites in the model, or a vector of indexes
%     for metabolites to exclude from this analysis (default []).
% isNames : logical
%     true if the supplied mets represent metabolite names (as opposed to
%     IDs). This is a way to delete metabolites in several compartments at
%     once without knowing the exact IDs. This only works if ignoreMets is a
%     cell array (default false).
% minNrFluxes : logical
%     solves the MILP problem of minimizing the number of fluxes instead of
%     the sum. Slower, but can be used if the sum gives too many fluxes
%     (default false).
% allowExcretion : logical
%     allow for excretion of all metabolites. Only used when direction is
%     'produce' (default true).
% params : struct
%     *obsolete option*.
% ignoreIntBounds : logical
%     true if internal bounds (including reversibility) should be ignored.
%     Exchange reactions are not affected. This can be used to find
%     unbalanced solutions which are not possible using the default
%     constraints (default false).
%
% Returns
% -------
% solution : double
%     flux vector for the solution.
% metabolite : double
%     the index of the metabolite(s) which was produced or consumed. If
%     possible only one metabolite is reported, but there are situations
%     where metabolites can only be exchanged in pairs (or more).
%
% Examples
% --------
%     [solution, metabolite] = findLeakMetabolite(model, 'produce');
%     [solution, metabolite] = findLeakMetabolite(model, 'consume', ignoreMets);
%
% Notes
% -----
% Works by forcing at least 1 unit of "any metabolite" to be produced or
% consumed, then minimising for the sum of fluxes. If more than one
% metabolite is produced/consumed, it tries to find a single metabolite
% that satisfies the constraint on its own.
%
% See Also
% --------
% makeSomething, consumeSomething

direction=char(direction);
if ~any(strcmp(direction,{'produce','consume'}))
    error('RAVEN:badInput','direction must be ''produce'' or ''consume''');
end

p=parseRAVENargs(varargin, {'ignoreMets',[]; ...
    'isNames',false; ...
    'minNrFluxes',false; ...
    'allowExcretion',true; ...
    'params',[]; ...
    'ignoreIntBounds',false});
ignoreMets=p.ignoreMets;
isNames=p.isNames;
minNrFluxes=p.minNrFluxes;
allowExcretion=p.allowExcretion;
params=p.params;
ignoreIntBounds=p.ignoreIntBounds;

if isempty(ignoreMets)
    ignoreMets=[];
elseif ~islogical(ignoreMets) && ~isnumeric(ignoreMets)
    ignoreMets=convertCharArray(ignoreMets);
end
if isempty(params)
    params.relGap=0.8;
end

if isNames==true && ~isempty(ignoreMets)
    if iscellstr(ignoreMets)==false %#ok<ISCLSTR>
        error('Must supply a cell array of strings if isNames=true');
    end
end

if isNames==false
    indexesToIgnore=getIndexes(model,ignoreMets,'mets');
else
    indexesToIgnore=[];
    for i=1:numel(ignoreMets)
        indexesToIgnore=[indexesToIgnore;find(strcmp(ignoreMets(i),model.metNames))]; %#ok<AGROW>
    end
end

if strcmp(direction,'produce')
    % Allow for excretion of all metabolites
    if allowExcretion==true
        model.b=[model.b(:,1) inf(numel(model.mets),1)];
    end
    exchCoeff    = -1;  % excretion reaction: metabolite leaves
    fakeRowCoeff =  1;  % excretion reactions produce the fake metabolite
    fakeExchCoeff= -1;  % excrete the fake metabolite
else  % 'consume'
    exchCoeff    = +1;  % uptake reaction: metabolite enters
    fakeRowCoeff = -1;  % uptake reactions consume the fake metabolite
    fakeExchCoeff= +1;  % uptake the fake metabolite
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

%Add exchange reactions for all metabolites (sign depends on direction)
model.S=[model.S speye(nMets)*exchCoeff];

%Add a fake metabolite coupled to all exchange reactions
model.S=[model.S;[sparse(1,nRxns) ones(1,nMets)*fakeRowCoeff]];

%Zero out exchange reactions for ignored metabolites (keeps indexing intact)
model.S(:,indexesToIgnore+nRxns)=0;

%Add an exchange reaction for the fake metabolite itself
model.S(size(model.S,1),size(model.S,2)+1)=fakeExchCoeff;
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
    %Find which exchange reactions are active; exclude the last fake one
    I=setdiff(find(sol.x(nRxns+1:end)>0.1),size(model.S,1));

    if any(I)
        %Try to narrow down to a single metabolite for interpretability
        oldS=model.S;
        foundSingle=false;
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
        if foundSingle==false
            model.S=oldS;
        end
        if minNrFluxes==true
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
