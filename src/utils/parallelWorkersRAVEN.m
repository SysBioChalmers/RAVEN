function nW = parallelWorkersRAVEN(runParallel)
% parallelWorkersRAVEN  Return the parfor worker count for a given setting.
%
% Returns a value for use with the parfor worker-count syntax:
%
%     nW = parallelWorkersRAVEN(runParallel);
%     parfor (i = 1:N, nW)
%         ...
%     end
%
% This avoids mutating global parallel-pool state: nW == 0 forces parfor
% to run as a plain for loop; nW == Inf lets the Parallel Computing Toolbox
% manage the pool size and auto-creation.
%
% Parameters
% ----------
% runParallel : logical
%     whether the caller should use parallel workers.
%
% Returns
% -------
% nW : double
%     0   to force sequential execution (parfor behaves as a for loop,
%         no pool required or created).
%     Inf to use available pool workers (PCT manages pool size and
%         auto-creation).

if ~runParallel
    nW = 0;
    return
end
if isempty(ver('parallel'))
    disp('Cannot find MATLAB Parallel Computing Toolbox, process is not parallelized.')
    nW = 0;
else
    nW = Inf;
end
end
