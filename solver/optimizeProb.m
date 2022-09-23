function res = optimizeProb(prob,params,verbose)
% optimizeProb
%   Optimize an LP or MILP formulated in cobra terms.
%
%   prob	cobra style LP/MILP problem struct to be optimised
%   params	solver specific parameters (optional)
%   verbose if true MILP progress is shown (opt, default true)
%
%   res		the output structure from the selected solver RAVENSOLVER
%   		(cobra style)

if nargin<2 || isempty(params)
    params=struct();
end
if nargin<3 || isempty(verbose)
    verbose = true;
end
%Set as global variable for speed improvement if optimizeProb is run many times
global RAVENSOLVER;
if isempty(RAVENSOLVER)
    if(~ispref('RAVEN','solver'))
        dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
    else
        RAVENSOLVER = getpref('RAVEN','solver');
    end
end
solver=RAVENSOLVER;



if ~all(lower(prob.vartype) == 'c')
    disp('MILP detected.');
    milp=true;
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
        solverparams.intTol = 10^-9; %min val for gurobi
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
    solverparams = structUpdate(solverparams,params);
    
    % Restructering problem according to gurobi format
    if isfield(prob, 'csense')
        prob.csense = renameparams(prob.csense, {'L','G','E'}, {'<','>','='});
        prob.sense = prob.csense;
        prob = rmfield(prob, {'csense'});
    end
    if isfield(prob, 'osense')
        prob.osense = renameparams(num2str(prob.osense), {'1','-1'}, {'min','max'});
        prob.modelsense = prob.osense;
        prob = rmfield(prob, {'osense'});
    end
    [prob.obj, prob.rhs] = deal(prob.c, prob.b);
    prob = rmfield(prob, {'c','b'});
    
    %Rename intTol to IntFeasTol
    if milp
        solverparams.IntFeasTol = solverparams.intTol;
        solverparams  = rmfield(solverparams, {'intTol'});
        prob.vtype = prob.vartype;
        prob  = rmfield(prob, {'vartype'});
    end
    
    resG = gurobi(prob,solverparams);
    
    try
        [res.full, res.obj, res.origStat] = deal(resG.x,  resG.objval, resG.status);
        if milp && strcmp(resG.status, 'TIME_LIMIT')
            % If res has the objval field, it succeeded, regardless of
            % time_limit status
            resG.status = 'OPTIMAL';
        end
        switch resG.status
            case 'OPTIMAL'
                res.stat = 1;
            case 'UNBOUNDED'
                res.stat = 2;
            otherwise
                res.stat = 0;
        end
        if ~milp
            [res.vbasis, res.cbasis] = deal(resG.vbasis, resG.cbasis);
		else
			res.mipgap = resG.mipgap; 
        end
    catch
        res.stat = 0;
        res.origStat = resG.status;  % useful information to have
    end
    %% Use GLPK using RAVEN-provided binary
    case 'glpk'
        solverparams.scale   = 128; % Auto scaling
        solverparams.tmlim   = defaultparams.timeLimit;
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
        solverparams.scale   = 1; % Auto scaling
        
        % Ensure that RAVEN glpk binary is used, return to original
        % directory afterwards
        [ravenDir,currDir]=findRAVENroot();
        cd(fullfile(ravenDir,'software','GLPKmex'))
        [xopt, fmin, errnum, ~] = glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
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
    otherwise
        error('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end
if res.stat>0
    res.full=res.full(1:size(prob.a,2));
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
