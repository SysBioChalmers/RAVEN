function res = optimizeProb(prob,params,verbose)
% optimizeProb
%   Optimize an LP or MILP formulated in cobra terms.
%
%   prob	cobra style LP/MILP problem struct to be optimised
%   params	solver specific parameters (optional)
%   verbose if true MILP progress is shown (opt, default true)
%
%   res		the output structure from the selected solver RAVENSOLVER
%   		(cobra style)

if nargin<2 || isempty(params)
    params=struct();
end
if nargin<3 || isempty(verbose)
    verbose = true;
end

%Set as global variable for speed improvement if optimizeProb is run many times
global RAVENSOLVER;
if isempty(RAVENSOLVER)
    if(~ispref('RAVEN','solver'))
        dispEM('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
    else
        RAVENSOLVER = getpref('RAVEN','solver');
    end
end
solver=RAVENSOLVER;

if ~all(lower(prob.vartype) == 'c')
    disp('MILP detected.');
    milp=true;
else
    milp=false;
end

%% Define default parameters, which will then be used to make solver-
% specific solverparams structures
defaultparams.feasTol        = 1e-9;
defaultparams.optTol         = 1e-9;
defaultparams.objTol         = 1e-6;
defaultparams.timeLimit      = 1000;
%defaultparams.iterationLimit = 1000;
defaultparams.intTol         = 1e-12;
defaultparams.relMipGapTol   = 1e-12;
defaultparams.absMipGapTol   = 1e-12;
if milp
    defaultparams.MIPGap     = 1e-12;
    defaultparams.Seed       = 1;
end

switch solver
    %% Use whatever solver is set by COBRA Toolbox changeCobraSolver
    case 'cobra'
        if milp
            cparams=struct('timeLimit',1e9,'printLevel',0,'intTol',1e-6,'relMipGapTol',1e-9);
            cparams=structUpdate(cparams,params);
            res=solveCobraMILP(prob,cparams);
        else
            res=solveCobraLP(prob);
        end
        if isfield(res,{'dual','rcost'})
            res.dual=res.dual;
            res.rcost=res.rcost;
        end

        %% Use Gurobi in a MATLAB environment
    case 'gurobi'
        if milp
            if verbose
                solverparams.OutputFlag = 1;
            else
                solverparams.OutputFlag = 0;
            end
            solverparams.IntFeasTol = 10^-9; %min val for gurobi
            solverparams.MIPGap = defaultparams.MIPGap;
            solverparams.Seed = defaultparams.Seed;
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
        if isfield(prob, 'csense')
            prob.sense = renameparams(prob.csense, {'L','G','E'}, {'<','>','='});
            prob = rmfield(prob, {'csense'});
        end
        if isfield(prob, 'osense')
            osense = prob.osense;
            prob.modelsense = renameparams(num2str(prob.osense), {'1','-1'}, {'min','max'});
            prob = rmfield(prob, {'osense'});
        end
        [prob.obj, prob.rhs, prob.vtype] = deal(prob.c, prob.b, prob.vartype);
        prob = rmfield(prob, {'c','b','vartype'});

        resG = gurobi(prob,solverparams);

        try
            % Name output fields the same as COBRA does
            res.full     = resG.x;
            res.obj      = resG.objval;
            res.origStat = resG.status;
            if isfield(resG,{'pi','rc'})
                res.dual     = -resG.pi*osense;
                res.rcost    = -resG.rc*osense;
            end
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
                res.vbasis = resG.vbasis;
                res.cbasis = resG.cbasis;
    		else
    			res.mipgap = resG.mipgap;
            end
        catch
            res.stat = 0;
            res.origStat = resG.status;  % useful information to have
        end
        %% Use GLPK using RAVEN-provided binary
    case 'glpk'
        solverparams.scale   = 1; % Auto scaling
        solverparams.tmlim   = defaultparams.timeLimit;
        solverparams.tolbnd  = defaultparams.feasTol;
        solverparams.toldj   = defaultparams.optTol;
        solverparams.tolint  = defaultparams.intTol;
        solverparams.tolobj  = defaultparams.objTol;
        solverparams.msglev  = 0; % Level of verbosity
        solverparams = structUpdate(solverparams,params);

        prob.csense = renameparams(prob.csense, {'L','G','E'}, {'U','L','S'});

        if milp
            solverparams.tmlim   = solverparams.tmlim*10;
            solverparams.msglev  = 1; % Level of verbosity
            disp('Issues have been observed when using GLPK for MILP solving. Be advised to carefully observe the results, or us another solver.')
        end

        % Ensure that RAVEN glpk binary is used, return to original
        % directory afterwards
        [ravenDir,currDir]=findRAVENroot();
        cd(fullfile(ravenDir,'software','GLPKmex'))
        [xopt, fmin, errnum, extra] = glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
        cd(currDir)

        switch errnum % 1 = undefined; 2 = feasible; 3 = infeasible; 4 = no feasible solution; 5 = optimal; 6 = no unbounded solution
            case 5
                res.stat = 1; % Optimal
            case 2
                res.stat = 2; % Feasible, but not optimal
            otherwise
                res.stat = 0;
        end
        res.origStat = errnum;
        res.full     = xopt;
        res.obj      = fmin;
        res.dual     = -extra.lambda*prob.osense;
        res.rcost    = -extra.redcosts*prob.osense;
        %% Use SoPlex
    case 'soplex'
        solverparams.tmlim   = 1;
        solverparams.save    = 1;
        solverparams = structUpdate(solverparams,params);
        prob.csense = renameparams(prob.csense, {'L','G','E'}, {'U','L','S'});

        if ismac
            try
                [~, zPath] = system("awk '/PATH=/{print}' ~/.zprofile");
                zPath = strip(regexprep(zPath,'^PATH="?',''));
                zPath = strsplit(regexprep(zPath,'"?',''),':');
                zPath(strcmp(zPath,'${PATH}')) = [];
                setenv('PATH', strcat(strjoin(zPath,':'), ':', getenv("PATH")));
            catch
            end
        end

        [ravenDir,currDir]=findRAVENroot();
        cd(fullfile(ravenDir,'software','GLPKmex'))
        glpk(prob.c, prob.A, prob.b, prob.lb, prob.ub, prob.csense, prob.vartype, prob.osense, solverparams);
        solverparams.tmlim   = defaultparams.timeLimit;
        [runCheck, cmdOut]  = system(['soplex --solvemode=2 -t' num2str(solverparams.tmlim) ' -x outpb.lp > result.out']);

        if runCheck ~= 0
            error('SoPlex did not run')
        end
        if verLessThan('matlab','9.9') %readlines introduced 2020b
            fid=fopen('result.out');
            line_raw=cell(1000000,1);
            i=1;
            while ~feof(fid)
                line_raw{i}=fgetl(fid);
                i=i+1;
            end
            line_raw(i:end)=[];
            line_raw=string(line_raw);
        else
            %line_raw=readlines('result.out','EmptyLineRule','skip','WhitespaceRule','trim');
            line_raw=readlines('result.out');
        end
        delete 'result.out'
        delete 'outpb.lp'
        cd(currDir)        
       
        if find(contains(line_raw,'problem is solved [optimal]'),1) > 0
            res.full = zeros(length(prob.c),1);
            res.stat = 1;
            obj_line = line_raw{contains(line_raw,'Objective value')};
            res.obj  = str2double(obj_line(strfind(obj_line,':')+2:end));
            a        = find(contains(line_raw,'Primal solution (name, value):'));
            z        = find(contains(line_raw,'All other variables are zero'));
            flux     = split(line_raw(a+1:z-1));
            rxns     = str2double(replace(flux(:,1),'x_',''));
            flux     = str2double(flux(:,2));
            res.full(rxns) = flux;
        elseif ~isempty(cmdOut)
            error(['SoPlex error: ' extractBefore(cmdOut,newline)])
        else
            statusLine = contains(line_raw,'SoPlex status');
            disp(line_raw(statusLine));
            res.stat = 0;
            res.obj  = [];
            res.full = [];
        end

    otherwise
        error('RAVEN solver not defined or unknown. Try using setRavenSolver(''solver'').');
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
