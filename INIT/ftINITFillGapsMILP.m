function [x,I,exitFlag]=ftINITFillGapsMILP(model, toMinimize, params, scores, verbose)
% ftINITFillGapsMILP
%   Returns the minimal set of fluxes that satisfy the model using
%   mixed integer linear programming. This is an optimized variant of the
%   old function "getMinNrFluxes" that is adapted to ftINIT.
%   It does not need to make an irrev model, which takes time. The problem 
%   also becomes smaller (fewer integers but larger matrix). Only tested with 
%   Gurobi.
%
%	model         a model structure
%   toMinimize    either a cell array of reaction IDs, a logical vector
%                 with the same number of elements as reactions in the model,
%                 of a vector of indexes for the reactions that should be
%                 minimized (optional, default model.rxns)
%   params        parameter structure as used by getMILPParams (optional)
%   scores        vector of weights for the reactions. Negative scores
%                 should not have flux. Positive scores are not possible in this
%                 implementation, and they are changed to max(scores(scores<0)).
%                 Must have the same dimension as toMinimize (find(toMinimize)
%                 if it is a logical vector) (optional, default -1 for all reactions)
%   verbose       if true, the MILP progression will be shown. 
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
% Usage: [x,I,exitFlag]=getMinNrFluxes(model, toMinimize, params, scores)

% glpk solver as implemented by COBRA does not work well for MILP.
global CBT_MILP_SOLVER
if strcmp(getpref('RAVEN','solver'),'cobra') && strcmp(CBT_MILP_SOLVER,'glpk')
    dispEM('The current solver is set to ''cobra'', while in COBRA the MILP solver has been set to ''glpk''. The COBRA implementation of glpk is not well suitable for solving MILPs. Please install the Gurobi or an alternative MILP solver.',true);
end

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

%The trick to make this possible using a reversible model is that for 
%reversible reactions, we do the following. Set the flux Vi == Vipos - Vineg.
%What happens then is that if the flux is positive, Vipos will have a nonzero value,
%and when the flux is negative, Vineg will have a positive value. In addition
%we can add an arbitrary constant to both Vipos and Vineg. For example, if the
%flux Vi is -1, Vineg can be 4 and Vipos 3. This is however not a problem, 
%because we can be sure that Vineg + Vipos >= abs(Vi). Since we are interested
%in forcing the ints to be on when there is a flux, it doesn't matter if we overestimate
%the flux! So, we can simply constrain the boolean Yi to Yi*Maxflux >= Vineg + Vipos.

%The matrix then becomes as this:
%         S       p n int b     var
%         SSSSSSSS        -
%         SSSSSSSS         -
%         SSSSSSSS          -
%         SSSSSSSS           -
%         SSSSSSSS            -
%         SSSSSSSS             -
%         -       1 -             
%               -  1 -
%           -          M         -
%             -         M         -
%                 - - M         -
%                  - -   M         -
% An example with 8 rxns and 6 metabolites. - means -1, M max flux, and S is the S matrix.
% 4 rxns are to be minimized(1,3,5,7) and 1,7 are reversible. The p and n 
% are the Vipos and Vineg variables (2 rxns of each). The ints are the Yi for
% the variables that are to be minimized (the rest of the rxns doesn't have any).
% The mets here are the constraints, so right under the S matrix, you have
% Vi == Vipos - Vineg for the reactions 1 and 7 while the two next rows represent
% the non-reversible rxns 3 and 5, where we simply say that yi*M >= Vi. The last 
% 2 rows are the reactions yi*M >= Vipos + Vineg. To the right, we first have a
% -I matrix for setting the b constraints, and under that we have rxns that are
% just variables (rxns) between 0 and Inf to complete the constraints mentioned above.
%Ex: yi*M >= Vipos + Vineg is impl. as yi*M - Vipos - Vineg - var == 0, 0 <= var <= Inf.
%
%All rows should be equal to zero, so we don't set the b vector in the problem
%The reactions should be constrained as follows
%S - as given in model.lb and model.ub
%pos and neg - between 0 and inf
%ints - between 0 and 1
%b - as stated in the model.b vector - if this has one column, that is lb and ub (fixed value), if two columns, that is lb and ub
%var - between zero and inf

[minLog, I]=ismember(model.rxns,toMinimize);
indexes=find(minLog);
revIndexes = find(minLog & (model.rev == 1));
irrevIndexes = find(minLog & (model.rev == 0));
revIndexesInInd = find(model.rev(indexes) == 1);
irrevIndexesInInd = find(model.rev(indexes) == 0);

%Add binary constraints in the following manner: -  Add one unique
%"metabolite" for each integer reaction as a substrate.
%   These metabolites can have net production
%-  Add reactions for the production of each of those metabolites. The
%   amount produced in one reaction unit must be larger than the largest
%   possible flux in the model (but not too large to avoid bad scaling)

maxFlux=1000;
%we build the total matrix as blocks: [S pos neg int b var]

%s block
intArray=speye(numel(model.rxns))*-1;
intArrayRev=intArray(revIndexes,:);
intArrayIrrev = intArray(irrevIndexes,:);
sBlock=[model.S;intArrayRev;intArrayIrrev;sparse(numel(revIndexes), numel(model.rxns))]; %the S matrix and what is below

%pos/neg blocks
revposorneg1 = sparse(numel(model.mets), numel(revIndexes));
revpos2 = speye(numel(revIndexes));
revneg2 = -revpos2;
revposneg3 = sparse(numel(irrevIndexes), numel(revIndexes));
revposneg4 = revneg2;
posBlock = [revposorneg1;revpos2;revposneg3;revposneg4];
negBlock = [revposorneg1;revneg2;revposneg3;revposneg4];

%int block
int1 = sparse(numel(model.mets), numel(indexes));
int2 = sparse(numel(revIndexes), numel(indexes));
tmpEye = speye(numel(indexes))*maxFlux;
int3 = tmpEye(irrevIndexesInInd,:);%we here select only the irrev indexes from the indexes
int4 = tmpEye(revIndexesInInd,:);
intBlock = [int1;int2;int3;int4];

%b block
b1 = speye(numel(model.mets))*-1;
b2 = sparse(numel(indexes) + numel(revIndexes), numel(model.mets));
bBlock = [b1;b2];

%var block
var1 = sparse(numel(model.mets), numel(indexes));
var2 = sparse(numel(revIndexes), numel(indexes));
tmpEye = speye(numel(indexes))*-1;
var3 = tmpEye(irrevIndexesInInd,:);
var4 = tmpEye(revIndexesInInd,:);
varBlock = [var1;var2;var3;var4];

prob.A = [sBlock posBlock negBlock intBlock bBlock varBlock];
prob.a = prob.A;%I think this is needed as well

prob.c=[zeros(numel(model.rxns),1);zeros(numel(revIndexes)*2,1); scores(:);zeros(numel(model.mets) + numel(indexes),1)]; %Minimize the sum of reaction scores for reactions that are on

%ub and lb, text copied from above
%S - as given in model.lb and model.ub
%pos and neg - between 0 and inf
%ints - between 0 and 1
%b - as stated in the model.b vector - if this has one column, that is lb and ub (fixed value), if two columns, that is lb and ub
%var - between zero and inf

if size(model.b,2)==2
    bub = model.b(:,2);
else
    bub = model.b(:,1);
end

prob.lb = [model.lb;zeros(numel(revIndexes)*2,1);zeros(numel(indexes),1);model.b(:,1);zeros(numel(indexes),1)];
prob.ub = [model.ub;inf(numel(revIndexes)*2,1);ones(numel(indexes),1);bub;inf(numel(indexes),1)];

prob.b=zeros(size(prob.a,1), 1);

%Use the output from the linear solution as starting point. Only the values
%for the integer variables will be used, but all are supplied.
intsIndexes = find(prob.c ~= 0);
%The start point is not important (solved quickly anyway), so just skip it.
%prob.sol.int.xx=zeros(numel(prob.c),1);
%prob.sol.int.xx(intsIndexes(sol.x(indexes)>10^-3))=1;%these doesn't work for gurobi anyways...
prob.x0=[];
prob.vartype=repmat('C', 1, size(prob.A,2));
prob.vartype(intsIndexes) = 'B';
prob.csense = repmat('E', 1, size(prob.A,1));
prob.osense=1; %minimize the objective

%prob=rmfield(prob,{'blx','bux','blc','buc'});
params.intTol = 10^-9; %experment with this value
params.TimeLimit = 300;
params.Seed = 26;%This is weird - although it says "optimal solution found", we can get different results with different
                 %values of the objective function, where one is more optimal than the other (pretty big difference...)
%params.CSClientLog = 3;%generates a warning in gurobi, but may be of interest for other solvers

% Optimize the problem
res = optimizeProb(prob,params,verbose);
isFeasible=checkSolution(res);

if ~isFeasible
    x=[];
    I=[];
    exitFlag=-1;
    return;
end

x=res.full(1:numel(model.rxns));%the fluxes
I=res.full(intsIndexes) > 10^-3;%The margin for integers in gurobi is 10^-5, not 10^-12 that was previously used! Use 10^-3 to have some margin!

tmp = res.full(intsIndexes);
sel = (tmp > 10^-12) & (tmp < 0.5);
if sum(sel) > 0
	%This may indicate that there is a problem with the tolerances in the solver
	disp(['ftINITFillGapsMILP: Some variables meant to be boolean in the MILP have intermediate values. Num vars: ' num2str(sum(sel))])
end

end
