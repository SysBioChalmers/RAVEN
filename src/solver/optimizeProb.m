function res = optimizeProb(prob,params)
% optimizeProb
%   Optimize an LP or MILP formulated in cobra terms.
%
%   prob	cobra style LP/MILP problem struct to be optimised
%   params	solver specific parameters (optional)
%
%   res		the output structure from the selected solver RAVENSOLVER
%   		(cobra style)
    
if nargin<2 || isempty(params)
    params=struct();
end
if(~ispref('RAVEN','solver'))
    dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
end
if ~all(lower(prob.vartype) == 'c')
    disp('MILP detected.');
    milp=true;
else
    milp=false;
end

% Define default parameters
defaultparams.feasTol        = 1e-6;
defaultparams.optTol         = 1e-6;
defaultparams.objTol         = 1e-6;
defaultparams.timeLimit      = 1000;
%defaultparams.iterationLimit = 1000;
defaultparams.intTol         = 1e-5;
defaultparams.relMipGapTol   = 1e-4;
defaultparams.absMipGapTol   = 1e-12;

solver=getpref('RAVEN','solver');
if endsWithOct(solver,'_octave')
    % Set parameters
    solverparams.scale   = 128; % Auto scaling
    solverparams.tmlim   = defaultparams.timeLimit;
    solverparams.tolbnd  = defaultparams.feasTol;
    solverparams.toldj   = defaultparams.optTol;
    solverparams.tolint  = defaultparams.intTol;
    solverparams.tolobj  = defaultparams.objTol;
    solverparams.msglev  = 0; % Level of verbosity    
    solverparams = structUpdate(solverparams,params);
    
    prob.csense = renameparams(prob.csense, {'L','G','E'}, {'U','L','S'});
    
    if strcmp(solver,'glpk_octave')
        % Run Octave internal glpk function
        if milp
            solverparams.tmlim   = solverparams.tmlim*10;
            solverparams.msglev  = 1; % Level of verbosity
            disp('Unreliable MILP solving via Octave GLPK. Be advised to use an alternative solver')
        end
        [xopt, fmin, errnum, extra] = glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
        switch errnum % Throw message for most common errors
            case 4
                error('GLPK error 4: Invalid bounds.')
            case 9
                error('GLPK error 9: Time limit exhausted.')
            case 5
                error('GLPK error 5: Solver failed.')
            otherwise
                if errnum~=0
                    error(['GLPK error ' num2str(errnum)])
                end
        end
        switch extra.status % 1 = undefined; 2 = feasible; 3 = infeasible; 4 = no feasible solution; 5 = optimal; 6 = no unbounded solution
            case 5
                res.stat = 1; % Optimal
            case 2
                res.stat = 2; % Feasible, but not optimal
            otherwise
                res.stat = 0;
        end
        res.origStat = extra.status;
        res.full     = xopt;
        res.obj      = fmin;
    else
        % Call alternative solver via command line
        % First run glpk for 1 ms, just to get the lp file that is used as input
        % for the alternative solver.
        solverparams.save = 1;
        solverparams.tmlim = 1;
        [~, ~, errnum, ~] = glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
        % Currently only gurobi implemented, can probably easily be expanded to
        % other solvers, as long as they can take lp file as input and sol file
        % as output.
        solFile = '_tmp.sol';
        logFile = '_tmp.log';
        if strcmp(solver,'gurobi_octave')
            % String to be parsed to command line
            clparams = [' OptimalityTol='   num2str(solverparams.toldj) ...
                        ' FeasibilityTol='  num2str(solverparams.tolbnd) ...
                        ' MIPGap='          num2str(defaultparams.relMipGapTol) ...
                        ' MIPGapAbs='       num2str(defaultparams.absMipGapTol) ...
                        ' IntFeasTol='      num2str(defaultparams.intTol) ...
                        ' ResultFile='      solFile ...
                        ' LogFile='         logFile ' .\outpb.lp'];
            if milp % Can take longer while no intermediate output is shown
                disp('MILP is being solved ...')
            end
            [~, cmdout]=system(['gurobi_cl' clparams]);
            delete('outpb.lp');
            
            errorCode = regexp(cmdout,'Error (\d+):','tokens','once');
            if ~isempty(errorCode)
                switch errorCode{1}
                    case 10009
                        error('Gurobi error 10009: Failed to obtain a valid license.');
                    case 10012
                        error('Gurobi error 10012: Failed to read the requested file.');
                    case 10013
                        error('Gurobi error 10013: Failed to write the requested file.');
                    otherwise
                        error(['Gurobi error ' num2str(errorCode{1})]);
                end
            end
            res=readLPsolution();
            % TODO: check for feasible (non-optimal), should be stat = 2
            if any(regexp(cmdout,'Optimal objective '))
                res.stat = 1;
            else
                res.stat = 0;
            end            
        end
    end
    
else
    % No octave solvers, use MATLAB APIs
    switch solver
        case 'gurobi_matlab'
            if milp
                solverparams.OutputFlag = 1;
            else
                solverparams.OutputFlag = 0;
            end
            solverparams.DisplayInterval= 1; % Level of verbosity
            solverparams.TimeLimit      = defaultparams.timeLimit;
            solverparams.FeasibilityTol = defaultparams.feasTol;
            solverparams.OptimalityTol  = defaultparams.optTol;
            solverparams.Presolve       = 2;
            solverparams = structUpdate(solverparams,params);
            
            % Restructering problem according to gurobi format
            prob.csense = renameparams(prob.csense, {'L','G','E'}, {'<','>','='});
            prob.osense = renameparams(num2str(prob.osense), {'1','-1'}, {'min','max'});
            
            [prob.obj, prob.rhs, prob.sense, prob.modelsense] = deal(prob.c, prob.b, prob.csense, prob.osense);
            prob = rmfield(prob, {'c','b','csense','osense'});
            
            resG = gurobi(prob,solverparams);
            
            try
                [res.full, res.obj, res.origStat] = deal(resG.x,  resG.objval, resG.status);
                if milp && strcmp(resG.status, 'TIME_LIMIT')
                    % If res has the objval field, it succeeded, regardless of
                    % time_limit status
                    resG.status = 'OPTIMAL';
                end 
                switch resG.status
                    case 'OPTIMAL'
                        res.stat = 1;
                    case 'UNBOUNDED'
                        res.stat = 2;
                    otherwise
                        res.stat = 0;        
                end
                if ~milp
                    [res.vbasis, res.cbasis] = deal(resG.vbasis, resG.cbasis);
                end
            catch
                res.stat = 0;
                res.origStat = resG.status;  % useful information to have
            end
        case 'cobra'
            if milp
                cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
                cparams=structUpdate(cparams,params);
                res=solveCobraMILP(prob,cparams);
            else
                res=solveCobraLP(prob);
            end
        otherwise
            error('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');     
    end
end

if res.stat>0
    res.full=res.full(1:size(prob.a,2));
end
end

function s_merged=structUpdate(s_old,s_new)
%Remove overlapping fields from first struct;
%Obtain all unique names of remaining fields;
%Merge both structs
s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
names = [fieldnames(s_merged); fieldnames(s_new)];
s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end

function paramlist = renameparams(paramlist,old,new)
if ~iscell(paramlist)
    wasNoCell = true;
    paramlist={paramlist};
else
    wasNoCell = false;
end
for i=1:numel(old)
    paramlist = regexprep(paramlist,old{i},new{i});
end
if wasNoCell
    paramlist=paramlist{1};
end
end
