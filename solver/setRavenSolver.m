function out = setRavenSolver(solver)
% setRavenSolver
%   Sets the solver for RAVEN and optionally saves a default value
%	to MATLAB startup.
%
%   solver		a string ('gurobi','mosek', 'cobra', ...)
%
%   Usage: setRavenSolver('gurobi')
%
%   Daniel Hermansson, 2016-10-10
%   Eduard Kerkhoven, 2016-10-22 - Use Matlab preferences for solver selection

if (~ischar(solver)) dispEM('Input should be a string.'); end
	
setpref('RAVEN','solver',solver)

end