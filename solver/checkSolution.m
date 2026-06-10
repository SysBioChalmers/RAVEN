function [isFeasible, isOptimal]=checkSolution(res)
% checkSolution  Check if a Mosek solution is feasible and optimal.
%
% Parameters
% ----------
% res : struct
%     the output structure from mosekopt.
%
% Returns
% -------
% isFeasible : logical
%     true if the solution is feasible.
% isOptimal : logical
%     true if the solution is optimal.
%
% Examples
% --------
%     [isFeasible, isOptimal] = checkSolution(res);
%
% Notes
% -----
% This function also throws an error if the license has expired.

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
