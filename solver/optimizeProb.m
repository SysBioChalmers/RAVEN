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

% Define default parameters, which will then be used to make solver-
% specific solverparams structures
defaultparams.feasTol        = 1e-6;
defaultparams.optTol         = 1e-6;
defaultparams.objTol         = 1e-6;
defaultparams.timeLimit      = 1000;
%defaultparams.iterationLimit = 1000;
defaultparams.intTol         = 1e-5;
defaultparams.relMipGapTol   = 1e-4;
defaultparams.absMipGapTol   = 1e-12;

solver=getpref('RAVEN','solver');
switch solver
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
        
    case 'cobra'
        if milp
            cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
            cparams=structUpdate(cparams,params);
            res=solveCobraMILP(prob,cparams);
        else
            res=solveCobraLP(prob);
        end
        
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
