function setRavenSolver(solver)
% setRavenSolver
%   Sets the solver for RAVEN and optionally saves a default value
%	to MATLAB startup.
%
%   solver		the name of the solver. The following options are available:
%                   -gurobi (all versions)
%                   -cobra. If this option is set, then RAVEN uses the current
%                   COBRA solver. Make sure that the one is set using
%                   changeCobraSolver
%
%   Usage: setRavenSolver(solver)
%
%   Simonas Marcisauskas, 2017-09-19
%

if (~ischar(solver)) dispEM('Input should be a string.'); end

if strcmpi(solver,'gurobi') || strcmpi(solver,'cobra') || strcmpi(solver,'none')
    setpref('RAVEN','solver',solver)
elseif strcmpi(solver,'mosek')
    error('MOSEK is no longer supported from RAVEN 2.3.0. Choose ''gurobi'' or ''cobra''.')
else
    EM='Such solver is not compatible with RAVEN. Available options are ''gurobi'' and ''cobra''';
    dispEM(EM);
end
