function removeRavenFromPath(verbose)
% removeRavenFromPath
%   This function removes RAVEN (and all subdirectories) from the MATLAB
%   path. With multiple RAVEN installations, this function only removes the
%   one that is highest up in the path. If afterwards any RAVEN
%   installations remain, a warning is given that the user might want to
%   rerun removeRavenFromPath to remove the additional RAVEN installations.
%
% Input:
%   verbose         logical, whether progress should be reported
% Usage: removeRavenFromPath()

% Get current RAVEN directory
ravenDir=findRAVENroot();

if contains(ravenDir,'MATLAB Add-Ons')
    error('RAVEN is installed as MATLAB Add-On. You should uninstall RAVEN from the Add-On Manager instead.')
end

% Get current paths
currPath = strsplit(path(),';');

% Select RAVEN (sub)dirs in path
ravenSubDirs = startsWith(currPath,ravenDir);
ravenSubDirs = currPath(ravenSubDirs);

% Remove from path
if ~isempty(ravenSubDirs)
    fprintf('Removing %s and subfolders from the MATLAB path\n',ravenDir)
    rmpath(ravenSubDirs{:});
    savepath;
end
try
    findRAVENroot();
    warning(['It appears that another RAVEN installation was present in the '...
          'MATLAB path. To also remove that one, you should rerun removeRavenFromPath().'])
catch
end