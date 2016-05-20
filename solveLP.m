function [solution hsSolOut prob]=solveLP(model,minFlux,params,hsSol)
% solveLP
%   Solves a linear programming problem
%
%   model         a model structure
%   minFlux       determines if a second optimization should be performed
%                 in order to get rid of loops in the flux distribution
%                 0: no such optimization is performed
%                 1: the sum of abs(fluxes) is minimized. This is the
%                 fastest way of getting rid of loops
%                 2: the square of fluxes is minimized. This tends to
%                 distribute fluxes across iso-enzymes, which results in a
%                 larger number of reactions being used
%                 3: the number of fluxes is minimized. This can result
%                 in the flux distributions that are the easiest to
%                 interpret. Note that this optimization can be very slow
%                 (opt, default 0)
%   params        parameter structure as used by getMILPParams. Only used when
%                 minFlux==3 (opt)
%   hsSol         hot-start solution for the LP solver. This can
%                 significantly speed up the process if many similar
%                 optimization problems are solved iteratively. Only used if
%                 minFlux is 0 or 1 (opt)
%
%   solution
%         f       objective value
%         x       primal (flux distribution)
%         stat    exit flag 
%                 1: the optimization terminated successfully
%                 0: the solution is feasible, but not necessarily optimal
%                -1: no feasible solution found
%                -2: solution found, but flux minimization failed
%         msg     string describing the status of the optimization
%   hsSolOut      solution to be used as hot-start solution (see the input
%                 parameters). Only used if minFlux is 0 or 1
%
%   Usage: [solution hsSolOut]=solveLP(model,minFlux,params,hsSol)
%
%   Rasmus Agren, 2013-07-16
%

if nargin<2
    minFlux=0;
end
if nargin<3
    params.relGap=0.4;
end
if nargin<4
    hsSol=[];
end

%Default return values
hsSolOut=[];
solution.x=[];
solution.f=[];
solution.stat=-1;

%Ignore the hot-start if the previous solution wasn't feasible
if isfield(hsSol,'prosta')
   if strfind(hsSol.prosta,'INFEASIBLE')
        hsSol=[];
   end
end
    
% Setup the problem to feed to MOSEK.
prob=[];
prob.c=model.c*-1;
prob.a=model.S;
prob.blc=model.b(:,1);
%If model.b has two column, then they are for lower/upper bound on the RHS
prob.buc=model.b(:,min(size(model.b,2),2));
prob.blx=model.lb;
prob.bux=model.ub;

%If hot-start should be used
if ~isempty(hsSol)
   prob.sol.bas=hsSol;
   params.MSK_IPAR_SIM_HOTSTART=1;
end

%Use MSK_OPTIMIZER_FREE_SIMPLEX. This should not be necessary, but I've
%noticed that the interior point solver is not as good at finding feasible
%solutions.

%%%%%%
params.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_FREE_SIMPLEX';
%[crap,res] = mosekopt('minimize echo(0)',prob,getMILPParams(params));
res = optimizeProb(prob,params);
%isFeasible=res.isFeasible; isOptimal=res.isOptimal;
%%%%%%

%Check if the problem was feasible and that the solution was optimal
[isFeasible isOptimal]=checkSolution(res);

%If the problem was infeasible using hot-start it is often possible to
%re-solve it without hot-start and get a feasible solution
if ~isFeasible && ~isempty(hsSol)
    prob.sol=rmfield(prob.sol,'bas');
    %[crap,res] = mosekopt('minimize echo(0)',prob,getMILPParams(params));
    res=optimizeProb(prob,params);
    [isFeasible isOptimal]=checkSolution(res);
end

%Return without solution if the problem was infeasible
if ~isFeasible
    solution.msg='The problem is infeasible';
    return;
end
if ~isOptimal
    solution.msg='The problem is feasible, but not necessarily optimal';
    solution.stat=0;
else
    %All is well
    solution.stat=1;
    solution.msg='Optimal solution found';
end

%Construct the output structure
if isfield(res.sol,'bas')
    solution.x=res.sol.bas.xx;
    if minFlux<=1;
        hsSolOut=res.sol.bas;
        if(isfield(res,'vbasis')) % gurobi uses vbasis and cbasis as hotstart
            hsSolOut.vbasis=res.vbasis;
            hsSolOut.cbasis=res.cbasis;
        end
    end
    solution.f=res.sol.bas.pobjval;
else
    %Interior-point. This is not used at the moment
    solution.x=res.sol.itr.xx;
    solution.f=res.sol.itr.pobjval;
end

%If either the square, the number, or the sum of fluxes should be minimized
%then the objective function value should be fixed before another
%optimization. It is not correct to fix the reactions which participate in
%the objective function to their values in solution.x, as there can be
%multiple solutions with the same objective function value. In addition, this
%approach could result in numerical issues when several fluxes are fixed.
%Instead a new "fake metabolite" is added to the problem. This metabolite
%is produced by each reaction with the stoichiometry that reaction has in
%the objective function. The equality constraint of that "fake metabolite"
%is then set to be at least as good as the objective function value.
if minFlux~=0
    model.S=[model.S;prob.c'];
    model.mets=[model.mets;'TEMP'];
    
    %If the constraint on the objective function value is exact there is a
    %larger risk of numerical errors. However for the quadratic fitting
    %intervals are not allowed
    if minFlux~=2
        if size(model.b,2)==1
            model.b=[model.b model.b];
        end
        if solution.f<0
            model.b=[model.b;[-inf solution.f*0.999999]];
        else
            model.b=[model.b;[-inf solution.f*1.000001]];
        end
    else
        model.b=[model.b;ones(1,size(model.b,2))*solution.f];
    end
    
    switch minFlux
        %The sum of fluxes should be minimized
        case 1
            %Convert the model to the irreversible format
            revRxns=find(model.rev);
            if ~isempty(revRxns)
               iModel=convertToIrrev(model);
            else
                iModel=model;
            end
            
            %Minimize all fluxes
            iModel.c(:)=-1;
            sol=solveLP(iModel);
            
            %Map back to reversible fluxes
            if sol.stat>=0
                solution.x=sol.x(1:numel(model.c));
                solution.x(revRxns)=solution.x(revRxns)-sol.x(numel(model.c)+1:end);
            else
                fprintf('WARNING: Could not solve the problem of minimizing the sum of fluxes. Uses output from original problem\n');
                solution.stat=-2;
            end
        %The square of fluxes should be minimized. This only works when
        %there is no interval on the mass balance constraints (model.b is a
        %vector)
        case 2
%         if size(model.b,2)==1
%             qsol=solveQP(model,model.rxns,zeros(numel(model.lb),1));
%             %There is a problem that the solver seldom converges totally in this
%             %kind of minimization. Print a warning but use the fluxes
%             if any(qsol.x)
%                 solution.x=qsol.x;
%                 if qsol.stat==-1
%                     fprintf('WARNING: The quadratic fitting did not converge\n');
%                 end
%             else
%                 fprintf('WARNING: Could not solve the problem of minimizing the square of fluxes. Uses output from linear program\n');
%             end
%         else
%         	fprintf('WARNING: Cannot minimize square of fluxes when size(model.b,2)==2. Uses output from linear program\n');
%         end
        fprintf('WARNING: Quadratic solver currently not working. Uses output from original problem\n');
        solution.stat=-2;
        %The number of fluxes should be minimized
        case 3
            [qx,I]=getMinNrFluxes(model,model.rxns,params);
            qx(~I)=0;
            if any(I)
                solution.x=qx;
            else
                fprintf('WARNING: Could not solve the problem of minimizing the number of fluxes. Uses output from linear program\n');
                solution.stat=-2;
            end
    end
end
end