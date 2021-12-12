function binEnd = binaryEnding()
% binaryEnding
%   Returns the appropriate default binary ending dependent on OS.

if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='.exe';
else
    EM='Unknown OS, exiting.';
    disp(EM);
    return
end
