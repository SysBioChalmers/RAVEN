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
%   Daniel Hermansson, 2016-04-01

if nargin<2
    params=[];
end

global RAVENSOLVER;

if(isempty('RAVENSOLVER'))
	dispEM('Raven solver not defined. Try using setRavenSolver("solver").');
end

milp=false;
if(isfield(prob,'ints')), disp('MILP detected.'); milp=true; end

switch RAVENSOLVER
	case 'gurobi' 
		gparams=struct('Presolve',2,'TimeLimit',1e9,'OutputFlag',0,'MIPGap',.01);
		%gparams=structUpdate(gparams,params);
		res = gurobi(mosekToGurobiProb(prob), gparams);

		res=gurobiToMosekRes(res,length(prob.c),milp);

	case 'cobra'
		if (milp)
			cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-4,'relMipGapTol',.01);
			%cparams=structUpdate(cparams,params);
			res=solveCobraMILP(mosekToCobraProb(prob),cparams);
		else
			% no hot start =/
			%cparams=struct('PrintLevel',0);
			res=solveCobraLP(mosekToCobraProb(prob));
		end
		res=cobraToMosekRes(res,length(prob.c),milp);

	case 'mosek'
		[crap,res] = mosekopt(['minimize echo(0)'],prob,getMILPParams(params));
	otherwise
		dispEM('Raven solver not defined or unknown. Try using setRavenSolver("solver").')
end

function s_merged=structUpdate(s_old,s_new)
	%// Remove overlapping fields from first struct%// Obtain all unique names of remaining fields,%// Merge both structs
	s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
	names = [fieldnames(s_merged); fieldnames(s_new)];
	s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end

end