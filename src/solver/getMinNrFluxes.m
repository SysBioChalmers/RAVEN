function [x,I,exitFlag]=getMinNrFluxes(model, varargin)
% getMinNrFluxes  Find the minimal set of fluxes that satisfy the model.
%
% Uses mixed integer linear programming to find the minimal set of fluxes
% that satisfy the model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% toMinimize : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of indexes
%     for the reactions that should be minimized (default model.rxns).
% params : struct
%     *obsolete option*.
% scores : double
%     vector of weights for the reactions. Negative scores should not have
%     flux. Positive scores are not possible in this implementation, and
%     they are changed to max(scores(scores<0)), or to 0 if scores has no
%     negative entry at all. Must have the same dimension as toMinimize
%     (find(toMinimize) if it is a logical vector) (default -1 for all
%     reactions).
% formulation : char
%     which MILP formulation to build (default 'irrev'):
%
%     - 'irrev' : convert the model to irreversible format and use one
%       binary per irreversible reaction, so a reversible reaction that is
%       minimized costs two binaries. The magnitude of the fluxes is
%       estimated with an LP first, which also provides a start solution.
%     - 'reversible' : keep the model reversible and use one binary per
%       reaction, at the cost of a larger constraint matrix. Fewer binaries
%       makes this the faster formulation on large problems, but it needs a
%       solver that handles MILPs well and no LP estimate of the flux
%       magnitudes is made, so an arbitrary large flux is assumed.
% verbose : logical
%     if true, the MILP progression will be shown (default false).
% resolveTies : logical
%     if true, pin the MILP's degenerate optimum to a canonical answer instead of
%     relying on Seed alone: hold the primary (scores) objective at its optimum,
%     then minimize the count of "on" reactions, then, among the sparsest, minimize
%     their summed reaction-id rank (lowest ids preferred). Two extra MILP solves;
%     best-effort, so a phase that does not converge within TimeLimit is dropped
%     with a warning rather than changing the result. Only implemented for
%     formulation='reversible' -- errors if combined with 'irrev'. Matches
%     raven-toolbox's resolve_ties (raven-gecko-parity#104) (default false).
%
% Returns
% -------
% x : double
%     the corresponding fluxes for the full model.
% I : double
%     the indexes of the reactions in toMinimize that were used in the
%     solution.
% exitFlag : double
%     exit status:
%
%     - 1 : optimal solution found
%     - -1 : no solution found, either because the problem is infeasible or
%       because the solver reached its time limit before finding one
%     - -2 : a solution was found but is not proven optimal, because the
%       solver stopped before reaching optimality. It is returned, but it
%       may use more fluxes than the minimum
%
% Examples
% --------
%     [x, I, exitFlag] = getMinNrFluxes(model, toMinimize, params, scores);
%
% Notes
% -----
% Uses 1000 mmol/gDW/h as an arbitary large flux. Could possibly cause
% problems if the fluxes in the model are larger than that.

p=parseRAVENargs(varargin, {'toMinimize',[]; 'params',[]; 'scores',[]; ...
    'formulation','irrev'; 'verbose',false; 'resolveTies',false});
toMinimize=p.toMinimize;
params=p.params;
scores=p.scores;
formulation=p.formulation;
verbose=p.verbose;
resolveTies=p.resolveTies;
if isempty(toMinimize)
    toMinimize=model.rxns;
elseif ~islogical(toMinimize) && ~isnumeric(toMinimize)
    toMinimize=convertCharArray(toMinimize);
else
    toMinimize=model.rxns(toMinimize);
end

if isempty(formulation)
    formulation='irrev';
end
formulation=lower(char(formulation));
if ~ismember(formulation,{'irrev','reversible'})
    EM='Valid options for formulation are "irrev" or "reversible"';
    error('RAVEN:badInput', '%s', EM);
end

%For passing parameters to the solver
if isempty(params)
    params=struct();
end

if isempty(scores)
    %It says that the default is -1, but that is to fit with other code
    scores=ones(numel(toMinimize),1)*1;
else
    if numel(scores)~=numel(toMinimize)
        EM='The number of scores must be the same as the number of reactions to minimize';
        error('RAVEN:badInput', '%s', EM);
    end

    %Change positive scores to have a small negative weight. This is a
    %temporary solution. If no score is negative there is nothing to clamp
    %to (max([]) errors here rather than filling in a value, since the
    %right-hand side is a computed empty array, not the literal [] MATLAB
    %treats as a delete), so clamp to 0 instead: every reaction is then
    %equally (un)weighted.
    negScores=scores<0;
    if any(negScores)
        scores(~negScores)=max(scores(negScores));
    else
        scores(~negScores)=0;
    end

    %It says that the default is -1, but that is to fit with other code
    scores=scores*-1;
end

if strcmp(formulation,'reversible')
    [x,I,exitFlag]=minNrFluxesReversible(model,toMinimize,scores,params,verbose,resolveTies);
else
    if resolveTies
        EM='resolveTies is only implemented for formulation=''reversible''';
        error('RAVEN:badInput', '%s', EM);
    end
    [x,I,exitFlag]=minNrFluxesIrrev(model,toMinimize,scores,params,verbose);
end
end

function [x,I,exitFlag]=minNrFluxesIrrev(model,toMinimize,scores,params,verbose)
%One binary per irreversible reaction. Reversible reactions to minimize are
%split by convertToIrrev and therefore carry two binaries each.

exitFlag=1;

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

%Return an empty solution if the non-constrained problem could not be solved
if isempty(sol.x)
    x=[];
    I=[];
    exitFlag=-1;
    return;
end

%Take the maximal times 5 to have a safe margin. If it is smaller than 1000,
%then use 1000 instead.
maxFlux=max(max(sol.x)*5,1000);

intArray=speye(numel(irrevModel.rxns))*-1;
intArray=intArray(indexes,:);
prob.a=[irrevModel.S;intArray];
a=[sparse(numel(irrevModel.mets),numel(indexes));speye(numel(indexes))*maxFlux];
prob.a=[prob.a a];
prob.ints.sub=numel(irrevModel.rxns)+1:numel(irrevModel.rxns)+numel(indexes);

prob.c=[zeros(numel(irrevModel.rxns),1);scores(:);zeros(size(prob.a,1),1)]; %Minimize the number of fluxes
prob.A=[prob.a -speye(size(prob.a,1))];
prob.blc=[irrevModel.b(:,1);zeros(numel(indexes),1)];
if size(irrevModel.b,2)==2
    prob.buc=[irrevModel.b(:,2);inf(numel(indexes),1)];
else
    prob.buc=[irrevModel.b(:,1);inf(numel(indexes),1)];
end
prob.blx=[irrevModel.lb;zeros(numel(indexes),1)];
prob.bux=[irrevModel.ub;ones(numel(indexes),1)];
prob.lb = [prob.blx; prob.blc];
prob.ub = [prob.bux; prob.buc];
prob.osense=1;
prob.csense=repmat('E', 1, size(prob.a,1),1);
prob.b=zeros(size(prob.a,1), 1);

%Use the output from the linear solution as starting point. Only the values
%for the integer variables will be used, but all are supplied.
prob.sol.int.xx=zeros(numel(prob.c),1);
prob.sol.int.xx(prob.ints.sub(sol.x(indexes)>10^-12))=1;
prob.x0=[];
prob.vartype=repmat('C', size(prob.A,2), 1);
prob.vartype(prob.ints.sub) = 'I'; % with .lb = 0 and .ub = 1, they are binary
% integers (glpk in octave only allows "continuous" or "", not "binary")
prob=rmfield(prob,{'blx','bux','blc','buc'});

% Optimize the problem
res = optimizeProb(prob,params,verbose);
[isFeasible, isOptimal]=checkSolution(res);

if ~isFeasible
    x=[];
    I=[];
    exitFlag=-1;
    if res.hitTimeLimit
        %Nothing is returned either way, but the caller reports -1 as "no
        %feasible solution exists", which is not what happened here.
        EM='Time limit reached before finding a solution. Try increasing the TimeLimit parameter.';
        warning('RAVEN:warning', '%s', EM);
    end
    return;
end
if res.hitTimeLimit || ~isOptimal
    %A feasible but suboptimal solution, i.e. the solver stopped before
    %reaching optimality. The solution is still returned, but must not be
    %reported as optimal: it may use more fluxes than the minimum.
    %Not every solver reports this in its status, gurobi for one presents a
    %MILP solution found before the limit as optimal, so res.hitTimeLimit is
    %consulted as well.
    exitFlag=-2;
end

xx=res.full(1:numel(irrevModel.rxns));
I=res.full(numel(xx)+1:end);

%Map back to original model from irrevModel
x=xx(1:numel(model.rxns));
if numel(irrevModel.rxns)>numel(model.rxns)
    x(model.rev~=0)=x(model.rev~=0)-xx(numel(model.rxns)+1:end);
end

I=ismember(toMinimize,strrep(irrevModel.rxns(indexes(I>10^-12)),'_REV',''));
end

function [x,I,exitFlag]=minNrFluxesReversible(model,toMinimize,scores,params,verbose,resolveTies)
%One binary per reaction, on the reversible model. Trades a larger
%constraint matrix for half the binaries on reversible reactions.

% glpk solver as implemented by COBRA does not work well for MILP.
global CBT_MILP_SOLVER
if strcmp(getpref('RAVEN','solver'),'cobra') && strcmp(CBT_MILP_SOLVER,'glpk')
    error('RAVEN:badInput', '%s', 'The current solver is set to ''cobra'', while in COBRA the MILP solver has been set to ''glpk''. The COBRA implementation of glpk is not well suitable for solving MILPs. Please install the Gurobi or an alternative MILP solver.');
end

exitFlag=1;

%The trick to make this possible using a reversible model is that for
%reversible reactions, we do the following. Set the flux Vi == Vipos - Vineg.
%What happens then is that if the flux is positive, Vipos will have a nonzero value,
%and when the flux is negative, Vineg will have a positive value. In addition
%we can add an arbitrary constant to both Vipos and Vineg. For example, if the
%flux Vi is -1, Vineg can be 4 and Vipos 3. This is however not a problem,
%because we can be sure that Vineg + Vipos >= abs(Vi). Since we are interested
%in forcing the ints to be on when there is a flux, it does not matter if we overestimate
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
% the variables that are to be minimized (the rest of the rxns does not have any).
% The mets here are the constraints, so right under the S matrix, you have
% Vi == Vipos - Vineg for the reactions 1 and 7 while the two next rows represent
% the non-reversible rxns 3 and 5, where we simply say that yi*M >= Vi. The last
% 2 rows are the reactions yi*M >= Vipos + Vineg. To the right, we first have a
% -I matrix for setting the b constraints, and under that we have rxns that are
% just variables (rxns) between 0 and Inf to complete the constraints mentioned above.
%Ex: yi*M >= Vipos + Vineg is impl. as yi*M - Vipos - Vineg - var == 0, 0 <= var <= Inf.
%
%All rows should be equal to zero, so we do not set the b vector in the problem
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

%The int columns are ordered as indexes, i.e. by position in the model,
%while the scores arrive ordered as toMinimize
scores=scores(I(indexes));

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

intsIndexes = find(prob.c ~= 0);
prob.x0=[];
prob.vartype=repmat('C', 1, size(prob.A,2));
prob.vartype(intsIndexes) = 'B';
prob.csense = repmat('E', 1, size(prob.A,1));
prob.osense=1; %minimize the objective

%Defaults that suit this formulation, without discarding what the caller set
if ~isfield(params,'intTol')
    params.intTol = 10^-9;
end
if ~isfield(params,'TimeLimit')
    params.TimeLimit = 300;
end
if ~isfield(params,'Seed')
    %Although the solver reports "optimal solution found", different seeds
    %can return solutions with noticeably different objective values, so a
    %fixed one is used to keep results reproducible
    params.Seed = 26;
end

% Optimize the problem
res = optimizeProb(prob,params,verbose);
[isFeasible, isOptimal]=checkSolution(res);

if ~isFeasible
    x=[];
    I=[];
    exitFlag=-1;
    if res.hitTimeLimit
        %Nothing is returned either way, but the caller reports -1 as "no
        %feasible solution exists", which is not what happened here.
        EM='Time limit reached before finding a solution. Try increasing the TimeLimit parameter.';
        warning('RAVEN:warning', '%s', EM);
    end
    return;
end
if res.hitTimeLimit || ~isOptimal
    %A feasible but suboptimal solution, i.e. the solver stopped before
    %reaching optimality. The solution is still returned, but must not be
    %reported as optimal: it may use more fluxes than the minimum.
    %Not every solver reports this in its status, gurobi for one presents a
    %MILP solution found before the limit as optimal, so res.hitTimeLimit is
    %consulted as well.
    exitFlag=-2;
end

if resolveTies
    %This MILP is highly degenerate: many reaction subsets reach the same
    %score optimum, and the solver returns an arbitrary one -- reproducible
    %for a fixed solver build + Seed, but fragile to a solver version bump
    %or an unrelated change upstream of the MILP (see raven-gecko-parity#104,
    %where this was demonstrated on genome-scale Human-GEM: the same Seed,
    %same problem, but a different MILP construction path flipped a tied
    %pair of reactions). This runs a lexicographic phase 2, holding the
    %score objective at its optimum with a floor constraint: first minimise
    %the count of "on" reactions (the sparsest optimum), then, among the
    %sparsest, minimise their summed reaction-id rank (prefer lower ids).
    %Ports raven-toolbox's _resolve_ties_fill (raven-toolbox#114).
    %
    %intCols identifies the int/binary variable columns directly from the
    %block layout (prob.c = [zeros(rxns); zeros(2*revIndexes); scores(:);
    %zeros(mets+indexes)]) rather than via find(prob.c~=0): a zero-scored
    %reaction would otherwise silently drop out of the column count and
    %misalign every later index into "indexes".
    intCols = numel(model.rxns) + 2*numel(revIndexes) + (1:numel(indexes));
    primaryObj = full(prob.c(intCols)' * res.full(intCols));
    tol = max(abs(primaryObj)*1e-7, 1e-7);

    floorRow = sparse(1,size(prob.A,2));
    floorRow(intCols) = prob.c(intCols)';
    tieProb = prob;
    tieProb.A = [prob.A; floorRow];
    tieProb.a = tieProb.A;
    tieProb.b = [prob.b; primaryObj+tol];
    tieProb.csense = [prob.csense 'L'];

    %An integer objective needs only a cheap absolute-gap proof, not a tiny
    %relative one -- RAVEN's own MIPGap default (1e-12) can be slow to prove
    %on a genome-scale problem where a gap below 1 already proves optimality
    %exactly (mirrors raven-toolbox's own MIPGap=0/MIPGapAbs=0.4 for the same
    %two phases).
    tieParams = params;
    tieParams.MIPGap = 0;
    tieParams.MIPGapAbs = 0.4;

    unproven = {};

    %Phase 2a: fewest "on" reactions among the score-optimal solutions.
    countObj = zeros(size(prob.c));
    countObj(intCols) = 1;
    tieProb.c = countObj;
    res2a = optimizeProb(tieProb,tieParams,verbose);
    [isFeasible2a, isOptimal2a] = checkSolution(res2a);
    if ~isOptimal2a
        unproven{end+1} = 'phase2a-parsimony';
    end
    if isFeasible2a
        kmin = full(countObj'*res2a.full);

        %Phase 2b: among the sparsest, prefer lower reaction ids -- a
        %deterministic tie-break independent of solver/seed.
        countRow = sparse(1,size(tieProb.A,2));
        countRow(intCols) = 1;
        tieProb2b = tieProb;
        tieProb2b.A = [tieProb.A; countRow];
        tieProb2b.a = tieProb2b.A;
        tieProb2b.b = [tieProb.b; kmin+0.5];
        tieProb2b.csense = [tieProb.csense 'L'];

        [~, sortOrder] = sort(model.rxns(indexes));
        ranks = zeros(numel(indexes),1);
        ranks(sortOrder) = (1:numel(indexes))';
        idObj = zeros(size(prob.c));
        idObj(intCols) = ranks;
        tieProb2b.c = idObj;
        res2b = optimizeProb(tieProb2b,tieParams,verbose);
        [isFeasible2b, isOptimal2b] = checkSolution(res2b);
        if ~isOptimal2b
            unproven{end+1} = 'phase2b-idrank';
        end
        if isFeasible2b
            %Adopt the tie-resolved solution; exitFlag still reflects
            %whether the *primary* (scores) objective was proven optimal,
            %same as without resolveTies.
            res = res2b;
        end
    end
    if ~isempty(unproven)
        EM = ['getMinNrFluxes tie resolution did not converge (' strjoin(unproven,', ') '): the selection among equally score-optimal solutions is itself an unproven incumbent, so resolveTies has reduced but not removed the run-to-run spread. Raise params.TimeLimit for a proven tie-break.'];
        warning('RAVEN:warning', '%s', EM);
    end
end

x=res.full(1:numel(model.rxns));%the fluxes

ints = res.full(intsIndexes);
%The margin for integers in gurobi is 10^-5; use 10^-3 to have some margin.
I=ismember(toMinimize, model.rxns(indexes(ints > 10^-3)));

sel = (ints > 10^-12) & (ints < 0.5);
if sum(sel) > 0
    %This may indicate that there is a problem with the tolerances in the solver
    disp(['getMinNrFluxes: Some variables meant to be boolean in the MILP have intermediate values. Num vars: ' num2str(sum(sel))])
end

end
