function [isFeasible isOptimal]=checkSolution(res)
%   Checks if a solution from Mosek is feasible and optimal
%
%   res             the output structure from mosekopt
%
%   isFeasible      true if the solution is feasible
%   isoptimal       true if the solution is optimal
%
%   This function also throws an error if the license has expired.
%
%   Usage: [isFeasible isOptimal]=checkSolution(res)
%
%   Rasmus Agren, 2013-07-05
%

if res.rcode==1001
	dispEM('The Mosek licence has expired');
end
if res.rcode==1008
	dispEM('The Mosek licence file is missing');
end
if res.rcode==1010
	dispEM('The Mosek licence used only supports small problems (up to 300 variables). Have you requested the correct licence?');
end
isFeasible=false;
isOptimal=false;
if isfield(res,'sol')
    if isfield(res.sol,'bas')
        %There are several types of infeasibilities, but I consider them
        %all to be the same
        if isempty(strfind(res.sol.bas.prosta,'INFEASIBLE'))
            isFeasible=true;
        end
        %There are several types of optimality, but I consider them all to
        %be the same
        if any(strfind(res.sol.bas.solsta,'OPTIMAL'))
            isOptimal=true;
        end
    else
        if isfield(res.sol,'int')
            %There are several types of infeasibilities, but I consider them
            %all to be the same
            if isempty(strfind(res.sol.int.prosta,'INFEASIBLE'))
                isFeasible=true;
            end
            %There are several types of optimality, but I consider them all to
            %be the same
            if any(strfind(res.sol.int.solsta,'OPTIMAL'))
                isOptimal=true;
            end
        else
            %This is when the interior point solver is used. That is currently
            %not the case
            return;
        end
    end
end
end