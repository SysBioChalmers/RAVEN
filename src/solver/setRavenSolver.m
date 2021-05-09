function setRavenSolver(solver)
% setRavenSolver
%   Sets the solver for RAVEN and optionally saves a default value
%	to MATLAB startup.
%
%   solver  string with solver name. The following options are available:
%           glpk    works currently only in Octave      
%           gurobi  works both in Octave and MATLAB
%           cobra   will only work in MATLAB, as COBRA is not Octave compatible.
%                   If this option is set, then RAVEN parses the problem to
%                   COBRA and uses the solver that is set by changeCobraSolver.
%
%   setRavenSolver considers whether Octave or MATLAB is used, affecting what
%   approach will be used by optimizeProb.
%
%   Usage: setRavenSolver(solver)

if (~ischar(solver)) dispEM('Input should be a string.'); end

switch solver
    case 'glpk'
        if isOctave
            setpref('RAVEN','solver','glpk_octave')
        else
            error('glpk is not yet set up for MATLAB. Use "cobra" as alternative.')
        end
    case 'gurobi'
        if isOctave
            setpref('RAVEN','solver','gurobi_octave')
        else
            setpref('RAVEN','solver','gurobi_matlab')
        end
    case 'cobra'
        if isOctave
            error('COBRA toolbox is not Octave compatible')
        else
            setpref('RAVEN','solver',solver)
        end
    case 'mosek'
        error('MOSEK support has been discontinued since RAVEN 2.3.0.')
otherwise
    error('Invalid solver defined')
end
end
