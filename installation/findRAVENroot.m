function [ravenPath, prevDir] = findRAVENroot()
% findRAVENroot
%   Finds the root of the RAVEN directory by searching for the path to
%   RAVEN.png. Can also record the current directory, in case a function will
%   use the ravenPath to navigate to a precise folder, and it should return to
%   the previous directory afterwards. See e.g. optimizeProb calling glpk.
%
%   The resolved path is cached for the rest of the MATLAB session: this is
%   called on every solver invocation, and re-reading the RAVEN.ravenPath
%   preference from disk every time is both wasteful and, under the sustained
%   call volume of a full test run, an intermittent source of failures from
%   MATLAB's own preference-file I/O. Run `clear findRAVENroot` after
%   changing that preference (checkRaven does this already) to pick
%   up the change without restarting MATLAB.

persistent cachedPath
prevDir = pwd();
if ~isempty(cachedPath)
    ravenPath = cachedPath;
    return;
end

% A stored preference is only trusted if it still points at a real RAVEN
% install; otherwise fall through to walking up from the currently
% executing copy of this file. Without this check, a stale preference
% left over from a different RAVEN checkout on the same machine silently
% resolves to that other copy's data, not the one actually running.
ravenPath = '';
if ispref('RAVEN','ravenPath')
    prefPath = getpref('RAVEN','ravenPath');
    if isfile(fullfile(prefPath,'RAVEN.png'))
        ravenPath = prefPath;
    end
end
if isempty(ravenPath)
    ST=dbstack('-completenames');
    ravenPath = ST(strcmp({ST.name},'findRAVENroot')).file;
    rootFound = 0;
    while rootFound == 0
        isRoot = isfile(fullfile(ravenPath,'RAVEN.png'));
        if isRoot
            rootFound = 1;
        else
            ravenPathOld = ravenPath;
            ravenPath = fileparts(ravenPath);
            if strcmp(ravenPathOld,ravenPath)
                error('Cannot find the RAVEN root directory. Make sure you have not removed the RAVEN.png file from your RAVEN installation.')
            end
        end
    end
end
cachedPath = ravenPath;
end
