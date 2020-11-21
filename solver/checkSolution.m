function [isFeasible, isOptimal]=checkSolution(res)
%   Checks if a solution from Mosek is feasible and optimal
%
%   res             the output structure from mosekopt
%
%   isFeasible      true if the solution is feasible
%   isoptimal       true if the solution is optimal
%
%   This function also throws an error if the license has expired.
%
%   Usage: [isFeasible, isOptimal]=checkSolution(res)

isFeasible=false;
isOptimal=false;

switch res.stat
    case 1
        isOptimal=true;
        isFeasible=true;
    case 2
        isOptimal=false;
        isFeasible=true;
    otherwise
        isOptimal=false;
        isFeasible=false;
end
end
