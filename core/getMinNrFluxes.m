function [x,I,exitFlag]=getMinNrFluxes(model, toMinimize, params,scores)
% getMinNrFluxes
%   Returns the minimal set of fluxes that satisfy the model using
%   mixed integer linear programming.
%
%	model         a model structure
%   toMinimize    either a cell array of reaction IDs, a logical vector
%                 with the same number of elements as reactions in the model,
%                 of a vector of indexes for the reactions that should be
%                 minimized (opt, default model.rxns)
%   params        parameter structure as used by getMILPParams (opt)
%   scores        vector of weights for the reactions. Negative scores
%                 should not have flux. Positive scores are not possible in this
%                 implementation, and they are changed to max(scores(scores<0)).
%                 Must have the same dimension as toMinimize (find(toMinimize)
%                 if it is a logical vector) (opt, default -1 for all reactions)
%
%   x             the corresponding fluxes for the full model
%   I             the indexes of the reactions in toMinimize that were used
%                 in the solution
%   exitFlag      1: optimal solution found
%                -1: no feasible solution found
%                -2: optimization time out
%
%   NOTE: Uses 1000 mmol/gDW/h as an arbitary large flux. Could possibly
%   cause problems if the fluxes in the model are larger than that.
%
%   Usage: [x,I,exitFlag]=getMinNrFluxes(model, toMinimize, params, scores)
%
%   Rasmus Agren, 2017-02-28
%

exitFlag=1;

if nargin<2
    toMinimize=model.rxns;
else
    if ~iscell(toMinimize)
        toMinimize=model.rxns(toMinimize);
    end
end

%For passing parameters to the solver
if nargin<3
    params=struct();
end

if nargin<4
    %It says that the default is -1, but that is to fit with other code
    scores=ones(numel(toMinimize),1)*1;
else
    if numel(scores)~=numel(toMinimize)
        EM='The number of scores must be the same as the number of reactions to minimize';
        dispEM(EM);
    end
    
    %Change positive scores to have a small negative weight. This is a
    %temporary solution.
    scores(scores>=0)=max(scores(scores<0));
    
    %It says that the default is -1, but that is to fit with other code
    scores=scores*-1;
end

%Check if the model is in irreversible format
if any(model.rev)
    %Convert the model to irreversible format
    irrevModel=convertToIrrev(model);
    
    %Find the indexes for the reactions in toMinimize
    [indexes, I]=ismember(strrep(irrevModel.rxns,'_REV',''),toMinimize);
else
    irrevModel=model;
    
    %Find the indexes for the reactions in toMinimize
    [indexes, I]=ismember(irrevModel.rxns,toMinimize);
end

indexes=find(indexes);
%Adjust scores to fit with reversible
scores=scores(I(indexes));

%Add binary constraints in the following manner: -  Add one unique
%"metabolite" for each integer reaction as a substrate.
%   These metabolites can have net production
%-  Add reactions for the production of each of those metabolites. The
%   amount produced in one reaction unit must be larger than the largest
%   possible flux in the model (but not too large to avoid bad scaling)

%Calculate a solution to the problem without any constraints. This is to
%get an estimate about the magnitude of fluxes in the model and to get a
%feasible start solution.
sol=solveLP(irrevModel,1);

%Return an empty solution if the non-constrained problem couldn't be solved
if isempty(sol.x)
    x=[];
    I=[];
    exitFlag=-1;
    return;
end

%Take the maximal times 5 to have a safe margin. If it's smaller than 1000,
%then use 1000 instead.
maxFlux=max(max(sol.x)*5,1000);

prob.c=[zeros(numel(irrevModel.rxns),1);scores(:)]; %Minimize the number of fluxes
prob.blc=[irrevModel.b(:,1);zeros(numel(indexes),1)];
if size(irrevModel.b,2)==2
    prob.buc=[irrevModel.b(:,2);inf(numel(indexes),1)];
else
    prob.buc=[irrevModel.b(:,1);inf(numel(indexes),1)];
end
prob.blx=[irrevModel.lb;zeros(numel(indexes),1)];
prob.bux=[irrevModel.ub;ones(numel(indexes),1)];

intArray=speye(numel(irrevModel.rxns))*-1;
intArray=intArray(indexes,:);
prob.a=[irrevModel.S;intArray];
a=[sparse(numel(irrevModel.mets),numel(indexes));speye(numel(indexes))*maxFlux];
prob.a=[prob.a a];
prob.ints.sub=numel(irrevModel.rxns)+1:numel(irrevModel.rxns)+numel(indexes);

%Use the output from the linear solution as starting point. Only the values
%for the integer variables will be used, but all are supplied.
prob.sol.int.xx=zeros(numel(prob.c),1);
prob.sol.int.xx(prob.ints.sub(sol.x(indexes)>10^-12))=1;

% Optimize the problem
res = optimizeProb(prob,params);
isFeasible=checkSolution(res);

if ~isFeasible
    x=[];
    I=[];
    exitFlag=-1;
    return;
end

xx=res.sol.int.xx(1:numel(irrevModel.rxns));
I=res.sol.int.xx(numel(xx)+1:end);

%Check if Mosek aborted because it reached the time limit
if strcmp('MSK_RES_TRM_MAX_TIME',res.rcode)
    exitFlag=-2;
end

%Map back to original model from irrevModel
x=xx(1:numel(model.rxns));
if numel(irrevModel.rxns)>numel(model.rxns)
    x(model.rev~=0)=x(model.rev~=0)-xx(numel(model.rxns)+1:end);
end

I=ismember(toMinimize,strrep(irrevModel.rxns(indexes(I>10^-12)),'_REV',''));
end
