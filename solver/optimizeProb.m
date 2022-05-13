function res = optimizeProb(prob,params)
% optimizeProb
%   Optimize an LP or MILP formulated in cobra terms.
%
%   prob	cobra style LP/MILP problem struct to be optimised
%   params	solver specific parameters (optional)
%
%   res		the output structure from the selected solver RAVENSOLVER
%   		(cobra style)

if nargin<2 || isempty(params)
    params=struct();
end
if(~ispref('RAVEN','solver'))
    dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end
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
defaultparams.objTol         = 1e-9;
defaultparams.timeLimit      = 1000;
%defaultparams.iterationLimit = 1000;
defaultparams.intTol         = 1e-12;
defaultparams.relMipGapTol   = 1e-12;
defaultparams.absMipGapTol   = 1e-12;

solver=getpref('RAVEN','solver');

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
        solverparams.OutputFlag = 1;
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
    prob.csense = renameparams(prob.csense, {'L','G','E'}, {'<','>','='});
    prob.osense = renameparams(num2str(prob.osense), {'1','-1'}, {'min','max'});
    
    [prob.obj, prob.rhs, prob.sense, prob.modelsense] = deal(prob.c, prob.b, prob.csense, prob.osense);
    prob = rmfield(prob, {'c','b','csense','osense'});
    
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
