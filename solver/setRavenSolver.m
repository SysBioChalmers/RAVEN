function setRavenSolver(solver)
% setRavenSolver
%   Sets the solver for RAVEN and optionally saves a default value
%	to MATLAB startup.
%
%   solver  string with solver name. The following options are available:
%           glpk    uses RAVEN-provided binaries
%           gurobi  requires a working Gurobi installation
%           cobra   requires a working COBRA Toolbox installation. RAVEN
%                   parses the problem to COBRA and uses the solver that is
%                   set by changeCobraSolver.
%
%   Usage: setRavenSolver(solver)
solver=char(solver);
switch solver
    case 'cobra'
        global CBT_LP_SOLVER;
        if isempty(CBT_LP_SOLVER)
            error('COBRA toolbox is not initialized, run ''initCobraToolbox()''')
        else
            setpref('RAVEN','solver','cobra');
        end
    case 'mosek'
        error('MOSEK support has been discontinued since RAVEN 2.3.0.')
    case {'glpk','gurobi','soplex'}
        setpref('RAVEN','solver',solver)
    otherwise
        error('Invalid solver defined')
end
global RAVENSOLVER
RAVENSOLVER = solver;
end
