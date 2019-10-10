function res = optimizeProb(prob,params)
% optimizeProb
%   Optimize an LP or MILP formulated in mosek terms.
%
%   prob	mosek style LP/MILP problem struct to be optimised
%   params	solver specific parameters (optional)
%
%   res		the output structure from the selected solver RAVENSOLVER
%   		(mosek style)
%
%	Eduard Kerkhoven, 2019-10-10
%

if nargin<2
    params=struct();
end


if(~ispref('RAVEN','solver'))
    dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end

milp=false;
if(isfield(prob,'ints')), disp('MILP detected.'); milp=true; end

solver=getpref('RAVEN','solver');
if strcmp(solver,'gurobi')
    gparams=struct('Presolve',2,'TimeLimit',1000,'OutputFlag',1,'MIPGap',1e-12,'Seed',0,'FeasibilityTol',1e-9,'OptimalityTol',1e-9);
    if (~milp) gparams.OutputFlag=0; end
    res = gurobi(cobraToGurobiProb(prob),gparams);
    res = gurobiToCobraRes(res);
elseif strcmp(solver,'cobra')
    if (milp)
        cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
        cparams=structUpdate(cparams,params);
        res=solveCobraMILP(prob,cparams);
    else
        res=solveCobraLP(prob);
    end
    %res=cobraToMosekRes(res,size(prob.a,2),milp);
else
    dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end
if res.stat>0
    res.full=res.full(1:size(prob.a,2));
end

    function s_merged=structUpdate(s_old,s_new)
        %// Remove overlapping fields from first struct%// Obtain all
        %unique names of remaining fields,%// Merge both structs
        s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
        names = [fieldnames(s_merged); fieldnames(s_new)];
        s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
    end

end
