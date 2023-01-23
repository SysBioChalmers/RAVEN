function removeRavenFromPath()
% removeRavenFromPath
%   This function removes all RAVEN directories and subdirectories from the
%   MATLAB path. This function only removes RAVEN from the MATLAB path, it
%   does not delete the RAVEN folder itself. If RAVEN was installed as
%   MATLAB Add-On, the user is prompted to instead uninstall RAVEN via the
%   Add-On Manager.
%
% Usage: removeRavenFromPath()

% Check for installation as Add-On
addList = matlab.addons.installedAddons;
if any(strcmp(addList.Name,'RAVEN Toolbox'))
    error(['RAVEN is installed as MATLAB Add-On. You should uninstall RAVEN via '...
        'the Add-On Manager instead. Once uninstalled, you may attempt to run '...
        'removeRavenFromPath again to remove any potential remaining RAVEN installations.'])
end

% Get current paths
currPath  = transpose(strsplit(path(),{';',':'}));
ravenPath = false(numel(currPath),1);
for i=1:numel(currPath)
    dirCont = ls(currPath{i});
    dirCont = cellstr(dirCont(3:end,:));
    if any(contains(dirCont,{'ravenCobraWrapper.m','makeFakeBlastStructure.m','setRavenSolver.m'})) % A few (hopefully) unique RAVEN functions
        ravenPath(i)=true;
    end
end

ravenPath = unique(regexprep(currPath(ravenPath),'(\\|\/)((external)|(struct_conversion)|(solver))',''));
addOnDir = contains(ravenPath,'MATLAB Add-Ons');
if any(addOnDir)
    warning(['RAVEN is installed as MATLAB Add-On at the following directory, but MATLAB '...
        'somehow does not recognize this installation. You are advised to manually delete '...
        'the directory to remove any RAVEN remnants: ' ravenPath{addOnDir}]);
end
ravenSubDirs = currPath(startsWith(currPath,ravenPath));
% Remove from path
if ~isempty(ravenSubDirs)
    rmpath(ravenSubDirs{:});
    savepath;
end
end
