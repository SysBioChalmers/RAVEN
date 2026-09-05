function res = optimizeProb(prob,varargin)
% optimizeProb  Optimize an LP or MILP formulated in COBRA terms.
%
% Parameters
% ----------
% prob : struct
%     COBRA-style LP/MILP problem struct to be optimised.
%
% Name-Value Arguments
% --------------------
% params : struct
%     solver-specific parameters. In addition to solver parameters, the
%     field "maxRatio" can be set to a number > 1 to improve numerical
%     conditioning: before solving (LP only) any column whose coefficients
%     span more than maxRatio is split via auxiliary metabolites/variables,
%     preserving the feasible region. See splitProbForConditioning. The
%     field is consumed here and not forwarded to the solver.
% verbose : logical
%     if true MILP progress is shown (default true).
%
% Returns
% -------
% res : struct
%     the output structure from the selected solver RAVENSOLVER (COBRA
%     style).
%
% See also
% --------
% splitProbForConditioning

p=parseRAVENargs(varargin, {'params',[]; 'verbose',true});
params=p.params;
if isempty(params)
    params=struct();
end
verbose=p.verbose;

%Set as global variable for speed improvement if optimizeProb is run many times
global RAVENSOLVER;
%CBT_MILP_SOLVER is a COBRA Toolbox global, set by changeCobraSolver, and
%has to be declared as such to be readable in the MILP check below
global CBT_MILP_SOLVER;
if isempty(RAVENSOLVER)
    if(~ispref('RAVEN','solver'))
        error('RAVEN:badInput', '%s', 'RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
    else
        RAVENSOLVER = getpref('RAVEN','solver');
    end
end
solver=RAVENSOLVER;
if ~all(lower(prob.vartype) == 'c')
    milp=true;
    errorText = 'glpk is not suitable for solving MILPs, ';
    switch solver
        case 'glpk'
            error([errorText 'select a different solver with setRavenSolver().'])
        case 'cobra'
            if strcmp(CBT_MILP_SOLVER,'glpk')
                error([errorText 'select a different solver with changeCobraSolver() or setRavenSolver().'])
            end
    end
else
    milp=false;
end

%% Define default parameters, which will then be used to make solver-
% specific solverparams structures
defaultparams.feasTol        = 1e-9;
defaultparams.optTol         = 1e-9;
defaultparams.objTol         = 1e-6;
defaultparams.timeLimit      = 1000;
%defaultparams.iterationLimit = 1000;
defaultparams.intTol         = 1e-12;
defaultparams.relMipGapTol   = 1e-12;
defaultparams.absMipGapTol   = 1e-12;
if milp
    defaultparams.MIPGap     = 1e-12;
    defaultparams.Seed       = 1;
end

%% Optionally improve numerical conditioning before solving. This is opt-in
% through params.maxRatio and only applied to LPs. Columns whose coefficients
% span more than maxRatio are split via auxiliary metabolites/variables, which
% preserves the feasible region but reduces the spread within each column. The
% added rows/columns are stripped from the result further down.
condApplied = false;
if isfield(params,'maxRatio')
    maxRatio = params.maxRatio;
    params   = rmfield(params,'maxRatio');
    if ~milp && isnumeric(maxRatio) && isscalar(maxRatio) && isfinite(maxRatio) && maxRatio>1
        [prob, condInfo] = splitProbForConditioning(prob, maxRatio);
        condApplied = condInfo.applied;
    end
end

res.obj=[];
switch solver
    %% Use whatever solver is set by COBRA Toolbox changeCobraSolver
    case 'cobra'
        if milp
            cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
            cparams=structUpdate(cparams,params);
            res=solveCobraMILP(prob,cparams);
        else
            res=solveCobraLP(prob);
        end
        %% Use Gurobi in a MATLAB environment
    case 'gurobi'
        if milp
            if verbose
                solverparams.OutputFlag = 1;
            else
                solverparams.OutputFlag = 0;
            end
            solverparams.IntFeasTol = 10^-9; %min val for gurobi
            solverparams.MIPGap = defaultparams.MIPGap;
            solverparams.Seed = defaultparams.Seed;
        else
            solverparams.OutputFlag = 0;
        end
        solverparams.DisplayInterval= 1; % Level of verbosity
        solverparams.TimeLimit      = defaultparams.timeLimit;
        solverparams.FeasibilityTol = defaultparams.feasTol;
        solverparams.OptimalityTol  = defaultparams.optTol;
        solverparams.Presolve       = 2;
        % getCurrentTask requires the Parallel Computing Toolbox; without it
        % the call itself errors, so guard it rather than requiring the
        % toolbox for every gurobi solve.
        if exist('getCurrentTask','file') && ~isempty(getCurrentTask) % If run in parallel, then one thread per gurobi
            solverparams.Threads=1;
        end
        solverparams = structUpdate(solverparams,params);

        % Restructering problem according to gurobi format
        if isfield(prob, 'csense')
            prob.sense = renameparams(prob.csense, {'L','G','E'}, {'<','>','='});
            prob = rmfield(prob, {'csense'});
        end
        % gurobi minimises unless told otherwise, so an absent osense means
        % minimisation; osense is also needed below to sign the duals
        osense = 1;
        if isfield(prob, 'osense')
            osense = prob.osense;
            if osense < 0
                prob.modelsense = 'max';
            else
                prob.modelsense = 'min';
            end
            prob = rmfield(prob, {'osense'});
        end
        [prob.obj, prob.rhs, prob.vtype] = deal(prob.c, prob.b, prob.vartype);
        prob = rmfield(prob, {'c','b','vartype'});

        resG = gurobi(prob,solverparams);

        % Name output fields the same as COBRA does. Every field below is
        % optional in the gurobi result: which ones are present depends on
        % the status and on the method gurobi chose (a barrier solve returns
        % no simplex basis, for instance). They are therefore read
        % individually rather than under one try/catch, which would have
        % reported a missing basis on an optimal solve as "infeasible".
        res.origStat = resG.status;
        if isfield(resG,'x')
            res.full = resG.x;
        end
        if isfield(resG,'objval')
            res.obj = resG.objval;
        end
        if isfield(resG,'pi') && isfield(resG,'rc')
            res.dual     = -resG.pi*osense;
            res.rcost    = -resG.rc*osense;
        end
        if milp && strcmp(resG.status, 'TIME_LIMIT') && isfield(resG,'x') && isfield(resG,'objval')
            % An incumbent was found before the time limit. It is usable but
            % not proven optimal, so it is reported as feasible rather than
            % optimal; callers that must not act on a suboptimal MILP
            % solution (ftINITFillGapsMILP, solveLP) key on that distinction.
            resG.status = 'SUBOPTIMAL';
        end
        switch resG.status
            case 'OPTIMAL'
                res.stat = 1;
            case {'UNBOUNDED','SUBOPTIMAL'}
                res.stat = 2;
            otherwise
                res.stat = 0;
        end
        if ~milp
            if isfield(resG,'vbasis')
                res.vbasis = resG.vbasis;
            end
            if isfield(resG,'cbasis')
                res.cbasis = resG.cbasis;
            end
        elseif isfield(resG,'mipgap')
            res.mipgap = resG.mipgap;
        end
        if res.stat>0 && ~isfield(res,'full')
            % A status that claims a solution without returning one
            res.stat = 0;
        end
        %% Use GLPK using RAVEN-provided binary
    case 'glpk'
        solverparams.scale   = 1; % Auto scaling
        %solverparams.tmlim   = defaultparams.timeLimit;
        solverparams.tolbnd  = defaultparams.feasTol;
        solverparams.toldj   = defaultparams.optTol;
        solverparams.tolint  = defaultparams.intTol;
        solverparams.tolobj  = defaultparams.objTol;
        solverparams.msglev  = 0; % Level of verbosity
        solverparams = structUpdate(solverparams,params);

        prob.csense = renameparams(prob.csense, {'L','G','E'}, {'U','L','S'});

        if milp
            solverparams.tmlim   = solverparams.tmlim*10;
            solverparams.msglev  = 1; % Level of verbosity
            disp('Issues have been observed when using GLPK for MILP solving. Be advised to carefully observe the results, or us another solver.')
        end

        % Ensure that RAVEN glpk binary is used, return to original
        % directory afterwards
        [ravenDir,currDir]=findRAVENroot();
        cd(fullfile(ravenDir,'software','GLPKmex'))
        [xopt, fmin, errnum, extra] = glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
        cd(currDir)

        switch errnum % 1 = undefined; 2 = feasible; 3 = infeasible; 4 = no feasible solution; 5 = optimal; 6 = no unbounded solution
            case 5
                res.stat = 1; % Optimal
            case 2
                res.stat = 2; % Feasible, but not optimal
            otherwise
                res.stat = 0;
        end
        res.origStat = errnum;
        res.full     = xopt;
        res.obj      = fmin;
        res.dual     = -extra.lambda*prob.osense;
        res.rcost    = -extra.redcosts*prob.osense;
        %% Use scip
    case {'soplex','scip'} % Old 'soplex' option also allowed
        % scip takes a row range [rl,ru] rather than a right-hand side plus
        % a sense, so csense has to be translated into that range. Passing
        % [-b,b] is only the same problem when every row is an equality with
        % b = 0, which is the case for solveLP but for none of the MILPs
        % built elsewhere in RAVEN.
        rl = full(prob.b(:));
        ru = rl;
        csense = prob.csense(:);
        rl(csense=='L') = -inf;   % A*x <= b
        ru(csense=='G') = inf;    % A*x >= b
        [xopt,fval,exitflag] = scip([], prob.c, prob.A, rl, ru, prob.lb, prob.ub, prob.vartype);

        %   [x,fval,exitflag,stats] = scip(H, f, A, rl, ru, lb, ub, xtype, sos, qc, nl, x0, opts)
        %
        %   Input arguments*:
        %       H - quadratic objective matrix (sparse, optional [NOT TRIL / TRIU])
        %       f - linear objective vector
        %       A - linear constraint matrix (sparse)
        %       rl - linear constraint lhs
        %       ru - linear constraint rhs
        %       lb - decision variable lower bounds
        %       ub - decision variable upper bounds
        %       xtype - string of variable integrality ("c" continuous, "i" integer, "b" binary)
        %       sos - SOS structure with fields type, index and weight (see below)
        %       qc - Quadratic Constraints structure with fields Q, l, qrl and qru (see below)
        %       nl - Nonlinear Objective and Constraints structure (see below)
        %       x0 - primal solution
        %       opts - solver options (see below)
        %
        %   Return arguments:
        %       x - solution vector
        %       fval - objective value at the solution
        %       exitflag - exit status (see below)
        %       stats - statistics structure
        %
        %   Option Fields (all optional, see also optiset for a list):
        %       solverOpts - specific SCIP options (list of pairs of parameter names and values)
        %       maxiter - maximum LP solver iterations
        %       maxnodes - maximum nodes to explore
        %       maxtime - maximum execution time [s]
        %       tolrfun - primal feasibility tolerance
        %       display - solver display level [0-5]
        %       probfile - write problem to given file
        %       presolvedfile - write presolved problem to file
        %
        %   Return Status, as returned by the MatlabSCIPInterface build that
        %   RAVEN provides and that setRavenSolver points at
        %   (github.com/scipopt/MatlabSCIPInterface). Verified against that
        %   binary with a bounded, an infeasible and an unbounded problem:
        %       1 - Problem Solved to Optimality
        %       2 - Problem is Infeasible
        %       3 - Problem is Unbounded
        %
        %   The OPTI Toolbox build of scip uses a different, longer set in
        %   which optimality is 11, infeasibility 12 and unboundedness 13,
        %   with the lower numbers meaning "a limit was reached". Those codes
        %   are accepted below as well, so either build reports optimality
        %   correctly; only OPTI's 1 ("User Interrupted") is ambiguous, and
        %   reading it as optimal is no worse than the previous mapping,
        %   which recognised no status from the shipped build at all.
        
        res.origStat = exitflag;
        res.full = xopt;
        res.obj  = fval;

        switch exitflag
            case {1, 11}                  % solved to optimality
                res.stat = 1;
            case {3, 13}                  % unbounded
                res.stat = 2;
            case {5, 6, 7, 8, 9, 10}      % OPTI: a limit was reached, so a
                res.stat = 2;             % solution may exist but is not proven optimal
            otherwise                     % 2/12 infeasible, 0 unknown, 4, 14
                res.stat = 0;
        end
    otherwise
        error('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end
if res.stat>0
    res.full=res.full(1:size(prob.a,2));
end
%Strip the auxiliary rows/columns that were added for conditioning, so that
%the result matches the dimensions of the original problem. res.full is
%already truncated to the original variables above (prob.a is left untouched).
if condApplied
    if isfield(res,'dual') && numel(res.dual)>condInfo.nOrigRows
        res.dual=res.dual(1:condInfo.nOrigRows);
    end
    if isfield(res,'rcost') && numel(res.rcost)>condInfo.nOrigCols
        res.rcost=res.rcost(1:condInfo.nOrigCols);
    end
end
end

function s_merged=structUpdate(s_old,s_new)
%Remove overlapping fields from first struct;
%Obtain all unique names of remaining fields;
%Merge both structs
s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
names = [fieldnames(s_merged); fieldnames(s_new)];
s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end

function paramlist = renameparams(paramlist,old,new)
if ~iscell(paramlist)
    wasNoCell = true;
    paramlist={paramlist};
else
    wasNoCell = false;
end
for i=1:numel(old)
    paramlist = regexprep(paramlist,old{i},new{i});
end
if wasNoCell
    paramlist=paramlist{1};
end
end
