function [newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(model,models,allowNetProduction,useModelConstraints,supressWarnings,rxnScores)
% fillGaps
%   Uses template model(s) to fill gaps in a model
%
%   model               a model structure that may contains gaps to be filled
%   models              a cell array of reference models or a model structure.
%                       The gaps will be filled using reactions from these models
%   allowNetProduction  true if net production of all metabolites is
%                       allowed. A reaction can be unable to carry flux because one of
%                       the reactants is unavailable or because one of the
%                       products can't be further processed. If this
%                       parameter is true, only the first type of
%                       unconnectivity is considered (optional, default false)
%   useModelConstraints true if the constraints specified in the model
%                       structure should be used. If false then reactions
%                       included from the template model(s) so that as many
%                       reactions as possible in model can carry flux
%                       (optional, default false)
%   supressWarnings     false if warnings should be displayed (optional, default
%                       false)
%   rxnScores           array with scores for each of the reactions in the
%                       reference model(s). If more than one model is supplied,
%                       then rxnScores should be a cell array of vectors.
%                       The solver will try to maximize the sum of the
%                       scores for the included reactions (optional, default
%                       is -1 for all reactions)
%
%   newConnected        cell array with the reactions that could be
%                       connected. This is not calulated if
%                       useModelConstraints is true
%   cannotConnect       cell array with reactions that could not be
%                       connected. This is not calculated if
%                       useModelConstraints is true
%   addedRxns           cell array with the reactions that were added from
%                       "models"
%   newModel            the model with reactions added to fill gaps
%   exitFlag            1: optimal solution found
%                      -1: no feasible solution found
%                      -2: optimization time out
%
%   This method works by merging the model to the reference model(s) and
%   checking which reactions can carry flux. All reactions that can't
%   carry flux are removed (cannotConnect).
%   If useModelConstraints is false it then solves the MILP problem of
%   minimizing the number of active reactions from the reference models
%   that are required to have flux in all the reactions in model. This
%   requires that the input model has exchange reactions present for the
%   nutrients that are needed for its metabolism. If useModelConstraints is
%   true then the problem is to include as few reactions as possible from
%   the reference models in order to satisfy the model constraints.
%   The intended use is that the user can attempt a general gap-filling using
%   useModelConstraint=false or a more targeted gap-filling by setting
%   constraints in the model structure and then use
%   useModelConstraints=true. Say that the user want to include reactions
%   so that all biomass components can be synthesized. He/she could then
%   define a biomass equation and set the lower bound to >0. Running this
%   function with useModelConstraints=true would then give the smallest set
%   of reactions that have to be included in order for the model to produce
%   biomass.
%
% Usage: [newConnected, cannotConnect, addedRxns, newModel, exitFlag]=...
%           fillGaps(model,models,allowNetProduction,useModelConstraints,...
%           supressWarnings,rxnScores,params)

%If the user only supplied a single template model
if ~iscell(models)
    models={models};
end

if nargin<3
    allowNetProduction=false;
end
if nargin<4
    useModelConstraints=false;
end
if nargin<5
    supressWarnings=false;
end
if nargin<6
    rxnScores=cell(numel(models),1);
    for i=1:numel(models)
        rxnScores{i}=ones(numel(models{i}.rxns),1)*-1;
    end
end
if isempty(rxnScores)
    rxnScores=cell(numel(models),1);
    for i=1:numel(models)
        rxnScores{i}=ones(numel(models{i}.rxns),1)*-1;
    end
end
if nargin<7
    params=struct();
end

if ~iscell(rxnScores)
    rxnScores={rxnScores};
end

models=models(:);
rxnScores=rxnScores(:);

%Check if the original model has an unconstrained field. If so, give a
%warning
if supressWarnings==false
    if isfield(model,'unconstrained')
        EM='This algorithm is meant to function on a model with exchange reactions for uptake and excretion of metabolites. The current model still has the "unconstrained" field';
        dispEM(EM,false);
    else
        if isempty(getExchangeRxns(model,'both'))
            fprintf('NOTE: This algorithm is meant to function on a model with exchange reactions for uptake and excretion of metabolites. The current model does not seem to contain any such reactions.\n');
        end
    end
end

%Simplify the template models to remove constrained rxns. At the same time,
%check that the id of the template models isn't the same as the model. That
%would cause an error further down
for i=1:numel(models)
    models{i}.rxnScores=rxnScores{i};
    models{i}=simplifyModel(models{i},false,false,true);
    if strcmpi(models{i}.id,model.id)
        EM='The reference model(s) cannot have the same id as the model';
        dispEM(EM);
    end
end

%This is a rather ugly solution to the issue that it's a bit tricky to keep
%track of which scores belong to which reactions. This requires that
%removeReactions and mergeModels are modified to check for the new field.
model.rxnScores=zeros(numel(model.rxns),1);

%First merge all models into one big one
allModels=mergeModels([{model};models],'metNames',true);

%Add that net production is ok
if allowNetProduction==true
    %A second column in model.b means that the b field is lower and upper
    %bound on the RHS.
    model.b=[model.b(:,1) inf(numel(model.mets),1)];
    allModels.b=[allModels.b(:,1) inf(numel(allModels.mets),1)];
end

if useModelConstraints==true
    newConnected={};
    cannotConnect={};
    
    %Check that the input model isn't solveable without any input
    sol=solveLP(model);
    if ~isempty(sol.f)
        addedRxns={};
        newModel=model;
        exitFlag=1;
        return;
    end
    
    %Then check that the merged model is solveable
    sol=solveLP(allModels);
    if isempty(sol.f)
        EM='There are no reactions in the template model(s) that can make the model constraints satisfied';
        dispEM(EM);
    end
    
    %Remove dead ends for speed reasons. This has to be done here and
    %duplicate below because there is otherwise a risk that a reaction
    %which is constrained to something relevant is removed
    allModels=simplifyModel(allModels,false,false,false,true,false,false,false,[],true);
    allModels.c(:)=0;
else
    %Remove dead ends for speed reasons
    allModels=simplifyModel(allModels,false,false,false,true,false,false,false,[],true);
    allModels.c(:)=0;
    
    %If model constraints shouldn't be used, then determine which reactions
    %to force to have flux
    
    %Get the reactions that can carry flux in the original model
    originalFlux=haveFlux(model,1);
    
    %For the ones that can't carry flux, see if they can do so in the
    %merged model
    toCheck=intersect(allModels.rxns(strcmp(allModels.rxnFrom,model.id)),model.rxns(~originalFlux));
    
    %Get the ones that still cannot carry flux
    I=haveFlux(allModels,1,toCheck);
    
    %Get the reactions that can't carry flux in the original model, but can
    %in the merged one
    K=toCheck(I);
    
    %This is a temporary thing to only look at the non-reversible rxns.
    %This is because all reversible rxns can have a flux in the
    %irreversible model format that is used by getMinNrFluxes
    [~, I]=ismember(K,model.rxns);
    K(model.rev(I)~=0)=[];
    
    %Constrain all reactions in the original model to have a flux
    allModels.lb(ismember(allModels.rxns,K))=0.1;
    
    %Return stuff
    newConnected=K;
    cannotConnect=setdiff(model.rxns(~originalFlux),newConnected);
end

%Then minimize for the number of fluxes used. The fixed rxns doesn't need
%to participate
templateRxns=find(~strcmp(allModels.rxnFrom,model.id));
[~, J, exitFlag]=getMinNrFluxes(allModels,templateRxns,params,allModels.rxnScores(templateRxns));

%Remove everything except for the added ones
I=true(numel(allModels.rxns),1);
I(templateRxns(J))=false;
addedModel=removeReactions(allModels,I,true,true,true);

newModel=mergeModels({model;addedModel},'metNames',true);
addedRxns=setdiff(newModel.rxns,model.rxns);
newModel=rmfield(newModel,'rxnScores');
end
