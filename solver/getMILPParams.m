function mosekParams=getMILPParams(params)
% getMILPParams
%   Returns a MOSEK parameter structure used when solving mixed-integer
%   linear programming problems
%
%	params          structure with one or more of the following fields
%       maxTime     maximal running time in minutes (opt, default 600)
%       relGap      maximal relative gap between integer and relaxed
%                   solution in order to be considered optimal (0.0-1.0)
%                   (opt, default 0.2)
%       printReport true if the results of the optimization should be
%                   displayed on the screen (opt, default false)
%
%   mosekParams     a parameter structure to be used with MOSEK
%
%   Usage: mosekParams=getMILPParams(params)
%
%   Rasmus Agren, 2014-05-08
%

if nargin<1
    params=[];
end

mosekParams=params;
mosekParams.MSK_DPAR_MIO_TOL_ABS_RELAX_INT=10^-9;
mosekParams.MSK_DPAR_MIO_TOL_REL_GAP=0.05;

%NOTE: These options were removed or renamed in Mosek 8. Should be investigated.
%mosekParams.MSK_DPAR_MIO_TOL_REL_RELAX_INT=10^-9;
%mosekParams.MSK_DPAR_MIO_TOL_X=10^-9;
mosekParams.MSK_DPAR_MIO_TOL_FEAS=10^-9;
mosekParams.MSK_DPAR_BASIS_TOL_S=10^-9;
mosekParams.MSK_DPAR_BASIS_TOL_X=10^-9;
mosekParams.MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP=10^-9;

%Get the mosek version. This is a bit problematic since the Mosek function
%for getting the version came in version 7.
if isfield(params,'presolve')
    mosekParams.MSK_DPAR_OPTIMIZER_MAX_TIME=params.presolve;
    mosekParams=rmfield(mosekParams,'presolve');
else
    if any(strfind(evalc('mosekopt info'),'MOSEK Version 7'))
        mosekParams.MSK_IPAR_PRESOLVE_USE=1;
    else
        mosekParams.MSK_DPAR_PRESOLVE_TOL_LIN_DEP=10^-9;
        %Turn off the presolve. This is because Mosek sometimes returns non-feasible
        %solutions because of problems with the presolver. Should check if version
        %is <6.0.0.147
        mosekParams.MSK_IPAR_PRESOLVE_USE=0;
    end
end

%Use a starting integer solution if supplied. This has no effect if no such
%solution is supplied
%mosekParams.MSK_IPAR_MIO_CONSTRUCT_SOL=1;

%10 hours as default upper time limit
mosekParams.MSK_DPAR_OPTIMIZER_MAX_TIME=10*60*60;

if isfield(params,'maxTime')
    mosekParams.MSK_DPAR_OPTIMIZER_MAX_TIME=params.maxTime*60;
    mosekParams=rmfield(mosekParams,'maxTime');
end
if isfield(params,'relGap')
    mosekParams.MSK_DPAR_MIO_TOL_REL_GAP=params.relGap;
    mosekParams=rmfield(mosekParams,'relGap');
end
if isfield(params,'printReport')
    mosekParams=rmfield(mosekParams,'printReport');
end
end
