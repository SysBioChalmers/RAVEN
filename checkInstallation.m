function [currVer, installType] = checkInstallation(developMode, checkBinaries)
% checkInstallation  DEPRECATED. Use checkRaven instead.
%
%   Forwards all arguments and outputs to checkRaven.
%
%   NOTE: this function is run before RAVEN has been added to the MATLAB
%   path, so it must not call any other RAVEN functions until checkRaven has
%   added RAVEN to the path itself. It cannot use the shared
%   utils/deprecationWarning for the same reason.
%
% Usage: [currVer, installType] = checkInstallation(developMode, checkBinaries)

warning('RAVEN:deprecated', ...
    'checkInstallation is deprecated and will be removed in the next major release. Use checkRaven instead.')
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
