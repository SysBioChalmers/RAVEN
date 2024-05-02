function [ps, oldPoolAutoCreate] = parallelPoolRAVEN(runParallel)
% handleParallelRAVEN
%   Called by RAVEN functions that support parallel processing, to confirm
%   whether the MATLAB Parallel Computing Toolbox is installed.
%   - The toolbox is installed, and runParallel == true, ==> a parallel
%     pool is started.
%   - The toolbox is installed, but runParallel == false, ==> the auto-
%     creation of a parallel pool is disabled, to prevent "parfor" in
%     the target function to start a pool anyway.
%   - The toolbox is not installed, and runParallel == true, ==> a warning
%     is displayed that parallel computer is not possible.
%   - The toolbox is not installed, and runParallel == false, ==> the
%     target runs as intended, as "parfor" will automatically run in serial
%     mode instead.
%
% Input:
%   runParallel         logical, whether the target function (which calls
%                       parallelPoolRAVEN) should be run in parallel (opt,
%                       default true)
%
% Output:
%   ps                  parallel settings structure that will be used by
%                       the target function    
%   oldPoolAutoCreate   logical, to reset the original ps.Pool.AutoCreate
%                       setting once the target function has finished
%
% Use: [ps, oldPoolAutoCreate] = parallelPoolRAVEN(runParallel)

if nargin<1 || isempty(runParallel)
    runParallel = true;
end

addonList = matlab.addons.installedAddons;
ps = []; oldPoolAutoCreate = [];
if ~any(strcmpi(addonList.Name,'Parallel Computing Toolbox'))
    if runParallel % User wants parallel, but will not be possible
        disp('Cannot find MATLAB Parallel Computing Toolbox, process is not parallelized.')
    end
else
    if ~runParallel % User has Parallel toolbox, but does not want pool to start.
        % If pool is already running, then it will use the max workers
        % anyway, but this is probably okay.
        ps = parallel.Settings;
        oldPoolAutoCreate = ps.Pool.AutoCreate;
        ps.Pool.AutoCreate = false;
    else
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('IdleTimeout',120);
        end
    end
end
end