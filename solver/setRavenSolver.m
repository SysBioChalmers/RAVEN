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
% Usage: setRavenSolver(solver)
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
    case {'glpk','gurobi'}
        setpref('RAVEN','solver',solver)
    case 'soplex'
        error('The ''soplex'' solver is replaced by ''scip''.')
    case 'scip'
        if ~ispc
            try
                scip; % The user might have installed SCIP manually
                setpref('RAVEN','solver','scip')
            catch
                error('SCIP not found. RAVEN only provides the precompiled SCIP MEX binary for Windows. Instructions on how to compile the SCIP MEX file are found at https://github.com/scipopt/MatlabSCIPInterface')
            end
        else
            ravenDir = findRAVENroot();
            if ~exist(fullfile(ravenDir,'software','scip','scip.mexw64'),'file')
                try
                    disp('Downloading and installing RAVEN-provided SCIP MEX binary...')
                    websave(fullfile(ravenDir,'software','scip_mex_win.zip'),'https://github.com/SysBioChalmers/RAVEN/releases/download/v2.8.6/scip_mex_win.zip'); % Should be updated to release with SCIP functionality
                catch
                    error('Unable to download SCIP MEX binary from the RAVEN GitHub page.')
                end
                unzip(fullfile(ravenDir,'software','scip_mex_win.zip'),fullfile(ravenDir,'software','scip'));
                delete(fullfile(ravenDir,'software','scip_mex_win.zip'));
            end
            try
                scip; % The pre-compiled MEX file might fail
                setpref('RAVEN','solver','scip')
            catch
                error('Unable to use the RAVEN-provided precompiled SCIP MEX binary. Instructions on how to compile the SCIP MEX file are found at https://github.com/scipopt/MatlabSCIPInterface')
            end
        end
    otherwise
        error('Invalid solver defined')
end
global RAVENSOLVER
RAVENSOLVER = solver;
end
