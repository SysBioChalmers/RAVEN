function [ravenPath, prevDir] = findRAVENroot()
% findRAVENroot
%   Finds the root of the RAVEN directory, first by  by searching for the path to
%   RAVEN.png. Can also record the current directory, in case a function will
%   use the ravenPath to navigate to a precise folder, and it should return to
%   the previous directory afterwards. See e.g. optimizeProb calling glpk.

ST=dbstack('-completenames');
prevDir = pwd();
% A stored preference is only trusted if it still points at a real RAVEN
% install. Without this check, a stale preference left over from a
% different RAVEN checkout on the same machine would silently resolve to
% that other copy's data, not the one actually running; erroring instead
% surfaces the stale preference itself, rather than silently falling back
% to a directory the caller did not ask for.
ravenPath = '';
if ispref('RAVEN','ravenPath')
    prefPath = getpref('RAVEN','ravenPath');
    if ~isfile(fullfile(prefPath,'RAVEN2.png'))
        error(['The RAVEN root preference points to ' strrep(prefPath,'\','/') ...
            ', which does not contain RAVEN2.png. Update it with ' ...
            'setpref(''RAVEN'',''ravenPath'',<path>), or clear it with ' ...
            'rmpref(''RAVEN'',''ravenPath'') to resolve the root from the ' ...
            'currently executing copy of RAVEN instead.'])
    end
    ravenPath = prefPath;
end
if isempty(ravenPath)
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