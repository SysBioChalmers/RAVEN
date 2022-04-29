%Part of ftINIT, optimized version of the previous function from RAVEN.
function [addedRxns, newModel, exitFlag]=fillGapsOpt(tModel, origModel, tRefModel,allowNetProduction,supressWarnings,rxnScores,params)
% fillGapsOpt
%   Variant of fillGaps specially adapted to speed up generation of ftINIT models
%
%   tModel              model that contains the task-specific rxns.
%   origModel           model without task-specific rxns.         
%   tRefModel           reference tModel - the full tModel, containing all rxns
%                       used for gap-filling of tModel + the task-specific rxns
%   allowNetProduction  true if net production of all metabolites is
%                       allowed. A reaction can be unable to carry flux because one of
%                       the reactants is unavailable or because one of the
%                       products can't be further processed. If this
%                       parameter is true, only the first type of
%                       unconnectivity is considered (opt, default false)
%   useModelConstraints true if the constraints specified in the tModel
%                       structure should be used. If false then reactions
%                       included from the template tModel(s) so that as many
%                       reactions as possible in tModel can carry flux
%                       (opt, default false)
%   supressWarnings     false if warnings should be displayed (opt, default
%                       false)
%   rxnScores           scores for each of the reactions in the
%                       reference tModel. 
%                       The solver will try to maximize the sum of the
%                       scores for the included reactions (opt, default
%                       is -1 for all reactions)
%   params              parameter structure as used by getMILPParams (opt)
%
%   addedRxns           the rxns added
%   newModel            the tModel with reactions added to fill gaps
%   exitFlag            1: optimal solution found
%                      -1: no feasible solution found
%                      -2: optimization time out
%
%   This method works by merging the tModel to the reference model and
%   checking which reactions can carry flux. All reactions that can't
%   carry flux are removed.
%   If useModelConstraints is false it then solves the MILP problem of
%   minimizing the number of active reactions from the reference model
%   that are required to have flux in all the reactions in model. This
%   requires that the input tModel has exchange reactions present for the
%   nutrients that are needed for its metabolism. If useModelConstraints is
%   true then the problem is to include as few reactions as possible from
%   the reference models in order to satisfy the tModel constraints.
%   The intended use is that the user can attempt a general gap-filling using
%   useModelConstraint=false or a more targeted gap-filling by setting
%   constraints in the model structure and then use
%   useModelConstraints=true. Say that the user want to include reactions
%   so that all biomass components can be synthesized. He/she could then
%   define a biomass equation and set the lower bound to >0. Running this
%   function with useModelConstraints=true would then give the smallest set
%   of reactions that have to be included in order for the tModel to produce
%   biomass.
%
%   Usage: [newModel, exitFlag]=...
%           fillGaps(tModel,models,allowNetProduction,...
%           supressWarnings,rxnScores,params)

if nargin<4
    allowNetProduction=false;
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

%Check if the tModel has an unconstrained field. If so, give a
%warning - we should check if this should be here later.
%We don't need to check this here, it is checked further up
%if supressWarnings==false
%    if isfield(tModel,'unconstrained')
%        EM='This algorithm is meant to function on a tModel with exchange reactions for uptake and excretion of metabolites. The current tModel still has the "unconstrained" field';
%        dispEM(EM,false);
%    else
%        if isempty(getExchangeRxns(tModel,'both'))
%            fprintf('NOTE: This algorithm is meant to function on a tModel with exchange reactions for uptake and excretion of metabolites. The current tModel does not seem to contain any such reactions.\n');
%        end
%    end
%end

%Simplify the template models to remove constrained rxns. At the same time,
%check that the id of the template models isn't the same as the tModel. That
%would cause an error further down
tRefModel.rxnScores=rxnScores;
tRefModel=simplifyModel(tRefModel,false,false,true);
%if strcmpi(tRefModel.id,tModel.id)
%    EM='The reference tModel cannot have the same id as the tModel';
%    dispEM(EM);
%end
%%%%tIrrevRefModel.rxnScores = rxnScoresIrrev;
%%%%tIrrevRefModel=simplifyModel(tIrrevRefModel,false,false,true);

%This is a rather ugly solution to the issue that it's a bit tricky to keep
%track of which scores belong to which reactions. This requires that
%removeReactions and mergeModels are modified to check for the new field.
tModel.rxnScores=zeros(numel(tModel.rxns),1);

%First merge all models into one big one
fullModel = tRefModel;%mergeModels([{tModel};models],'metNames',true);

%temp for testing
%global m1;
%m1 = fullModel;

%Add that net production is ok
if allowNetProduction==true
    %A second column in tModel.b means that the b field is lower and upper
    %bound on the RHS.
    tModel.b=[tModel.b(:,1) inf(numel(tModel.mets),1)];
    fullModel.b=[fullModel.b(:,1) inf(numel(fullModel.mets),1)];
%%%%    tIrrevRefModel.b=[tIrrevRefModel.b(:,1) inf(numel(tIrrevRefModel.mets),1)];
end

% if useModelConstraints == true % we only use this case, the previously exising else was removed

%Check that the input tModel isn't solveable without any input
%We have already checked this in fitTasksOpt, no need to check it again
%sol=solveLP(tModel);
%if ~isempty(sol.f)
%    addedRxns={};
%    newModel=tModel;
%    exitFlag=1;
%    return;
%end

%Then check that the merged tModel is solveable
%we just assume it is - otherwise there is either an error in the original model or in the tasks.
%If there is problems with this, instead run the ordinary checkTasks and investigate this before using these with tINIT.
%sol=solveLP(fullModel);
%if isempty(sol.f)
%    EM='There are no reactions in the template tModel(s) that can make the tModel constraints satisfied';
%    dispEM(EM);
%end

%Remove dead ends for speed reasons. This has to be done here and
%duplicate below because there is otherwise a risk that a reaction
%which is constrained to something relevant is removed
%It is probably difficult to get rid of this.
%fullModel = simplifyModel(fullModel,false,false,false,true,false,false,false,[],true);
fullModel.c(:)=0;
%%%%tIrrevRefModel.c(:)=0;
%end useModelConstraints 


%Then minimize for the number of fluxes used. The fixed rxns doesn't need
%to participate
templateRxns = find(~ismember(fullModel.rxns, tModel.rxns)); %Check if this is slow, in that case keep track of this in fitTasksOpt and send it in
%templateRxns=find(~strcmp(fullModel.rxnFrom,tModel.id));

%%%[~, J, exitFlag]=getMinNrFluxesOpt2(fullModel,tIrrevRefModel,templateRxns,params,fullModel.rxnScores(templateRxns));%only the scores from the template rxns are used, so the others doesn't matter
%%%[~, J, exitFlag]=getMinNrFluxes(fullModel,templateRxns,params,fullModel.rxnScores(templateRxns));%only the scores from the template rxns are used, so the others doesn't matter
[~, J, exitFlag]=getMinNrFluxesOpt(fullModel,templateRxns,params,fullModel.rxnScores(templateRxns));%only the scores from the template rxns are used, so the others doesn't matter

%Remove everything except for the added ones
%%I=true(numel(fullModel.rxns),1);
%%I(templateRxns(J))=false;
%%addedModel=removeReactions(fullModel,I,true,true,true); %check if this is slow - no, quick
%I=true(numel(fullModel.rxns),1);
%I(templateRxns(J))=false;
addedRxns = fullModel.rxns(templateRxns(J));
rxnsToAdd.rxns = addedRxns;
rxnsToAdd.equations = constructEquations(fullModel, addedRxns);
rxnsToAdd.ub = fullModel.ub(templateRxns(J));
rxnsToAdd.lb = fullModel.lb(templateRxns(J));
newModel = addRxns(origModel,rxnsToAdd, 3, [], true); %we ignore gene rules etc here, they are not needed - we regenerate the model in the end by subtracting rxns from the full model.
%ismember(rxnsToAdd.rxns, tModel.rxns)

%rxnsToRem = setdiff(tRefModel.rxns, [tModel.rxns;tRefModel.rxns(templateRxns(J))]);
%newModel = removeReactions(origRefModel, rxnsToRem, true);%remove rxns and mets, skip genes I think for now
%addedRxns = tRefModel.rxns(templateRxns(J));
%%addedRxns = addedModel.rxns;
%%newModel=mergeModels({origModel;addedModel},'metNames',true);
%addedRxns=setdiff(newModel.rxns,origModel.rxns);
%%newModel=rmfield(newModel,'rxnScores');
end
