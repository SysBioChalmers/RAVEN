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
%	Eduard Kerkhoven, 2016-10-22 - Use Matlab preferences for solver selection
%

if nargin<2 || isempty(params)
    params=struct();
end


if(~ispref('RAVEN','solver'))
    dispEM(['Raven solver not defined or unknown. Try using setRavenSolver(',char(39),'solver',char(39),').']);
end

milp=false;
if(isfield(prob,'ints')), disp('MILP detected.'); milp=true; end

solver=getpref('RAVEN','solver');
if strcmp(solver,'gurobi')
    gparams=struct('Presolve',2,'TimeLimit',1000,'OutputFlag',1,'MIPGap',1e-12,'Seed',0,'FeasibilityTol',1e-9,'OptimalityTol',1e-9);
    if (~milp)
        gparams.OutputFlag=0;
    end
    gparams=structUpdate(gparams,params);
    % remove some MOSEK-specific fields that generate warnings with Gurobi
    gparams=rmfield(gparams, intersect(fieldnames(gparams),{'MSK_IPAR_OPTIMIZER';'relGap'}));
    res = gurobi(mosekToGurobiProb(prob), gparams);
    res=gurobiToMosekRes(res,length(prob.c),milp);
elseif strcmp(solver,'cobra')
    if (milp)
        cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
        cparams=structUpdate(cparams,params);
        res=solveCobraMILP(mosekToCobraProb(prob),cparams);
    else
        % no hot start =/
        cparams=struct('PrintLevel',1);
        res=solveCobraLP(mosekToCobraProb(prob));
    end
    res=cobraToMosekRes(res,length(prob.c),milp);
    
elseif strcmp(solver,'mosek')
    if (milp)
        echo=0;
        if isfield(params,'printReport') && params.printReport==true
            echo='3';
        end
        [~,res] = mosekopt(['minimize echo(' echo ')'],prob,getMILPParams(params));
    else
        [~,res] = mosekopt('minimize echo(0)',prob);
    end
else
    dispEM(['Raven solver not defined or unknown. Try using setRavenSolver(',char(39),'solver',char(39),')']);
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
