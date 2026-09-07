function [currVer, installType] = checkInstallation(developMode, checkBinaries)
% checkInstallation
%   Deprecated, use checkRaven instead. Forwards all arguments and outputs.
%
%   NOTE: this function is run before RAVEN has been added to the MATLAB
%   path, so it must not call any other RAVEN functions until checkRaven has
%   added RAVEN to the path itself.
%
% Usage: [currVer, installType] = checkInstallation(developMode, checkBinaries)

warning('checkInstallation is deprecated and will be removed in a future release, use checkRaven instead.')
if nargin < 1
    developMode = false;
end
if nargin < 2
    checkBinaries = true;
end
if nargout > 0
    [currVer, installType] = checkRaven(developMode, checkBinaries);
else
    checkRaven(developMode, checkBinaries);
end
end
