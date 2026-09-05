function deprecationWarning(oldFunction, replacement)
% deprecationWarning  Warn once per session that a function is deprecated.
%
% Used by the wrappers in the deprecated/ folder. Each function warns the
% first time it is called in a MATLAB session and stays quiet afterwards,
% so a loop over thousands of reactions does not bury the console. Call
% "clear functions" to reset.
%
% Parameters
% ----------
% oldFunction : char
%     name of the deprecated function.
% replacement : char
%     what to call instead, written as it would be typed, e.g.
%     'findLeakMetabolite(model,''produce'',...)'.
%
% Examples
% --------
%     deprecationWarning('makeSomething','findLeakMetabolite(model,''produce'',...)');

persistent alreadyWarned

if isempty(alreadyWarned)
    alreadyWarned=containers.Map('KeyType','char','ValueType','logical');
end

oldFunction=char(oldFunction);
if isKey(alreadyWarned,oldFunction)
    return
end
alreadyWarned(oldFunction)=true;

warning('RAVEN:deprecated', ...
    ['%s is deprecated and will be removed in the next major release. ' ...
     'Use %s instead.'], oldFunction, char(replacement));
end
