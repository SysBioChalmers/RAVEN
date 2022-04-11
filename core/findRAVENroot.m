function [ravenPath, prevDir] = findRAVENroot()
% findRAVENroot
%   Finds the root of the RAVEN directory, by searching for the path to
%   RAVEN2.png. Can also record the current directory, in case a function will
%   use the ravenPath to navigate to a precise folder, and it should return to
%   the previous directory afterwards. See e.g. optimizeProb calling glpk.

ST=dbstack('-completenames');
prevDir = pwd();
if length(ST)>1
    ravenPath=ST(2).file; % In case findRAVENroot is run via another function
else
    ravenPath=ST(1).file;
end
rootFound = 0;
while rootFound == 0
    isRoot = exist(fullfile(ravenPath,'RAVEN2.png'),'file');
    if isRoot == 2
        rootFound = 1;
    else
        ravenPathOld = ravenPath;
        ravenPath = fileparts(ravenPath);
        if strcmp(ravenPathOld,ravenPath)
            error('Cannot find the RAVEN root directory. Make sure you have not removed the RAVEN2.png file from your RAVEN installation.')
        end
    end
end