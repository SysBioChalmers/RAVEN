function removeRavenFromPath(verbose)
% removeRavenFromPath
%   This function removes RAVEN (and all subdirectories) from the MATLAB
%   path. With multiple RAVEN installations, this function only removes the
%   one that is highest up in the path. If afterwards any RAVEN
%   installations remain, a warning is given that the user might want to
%   rerun removeRavenFromPath to remove the additional RAVEN installations.
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
currPath  = transpose(strsplit(path(),';'));
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
