function legacyMethodNotice(oldFunction, preferred)
% legacyMethodNotice  Point a caller at the recommended method, once.
%
% For functions that remain supported but are no longer the recommended way
% to do the job. Unlike deprecationWarning, this promises no removal: the
% function keeps working and keeps being maintained. It is a signpost, not a
% countdown.
%
% The notice is issued the first time the function is called in a MATLAB
% session and stays quiet afterwards, so a script that calls it in a loop
% does not bury the console. Call "clear functions" to reset. A caller that
% wants the legacy method deliberately can silence it with
% warning('off','RAVEN:legacyMethod').
%
% Parameters
% ----------
% oldFunction : char
%     name of the function being called.
% preferred : char
%     what to use instead for new work, written as it would be typed.
%
% Examples
% --------
%     legacyMethodNotice('getINITModel','ftINIT');
%
% See Also
% --------
% deprecationWarning

persistent alreadyNoticed

if isempty(alreadyNoticed)
    alreadyNoticed=containers.Map('KeyType','char','ValueType','logical');
end

oldFunction=char(oldFunction);
if isKey(alreadyNoticed,oldFunction)
    return
end
alreadyNoticed(oldFunction)=true;

warning('RAVEN:legacyMethod', ...
    ['%s is the legacy implementation and remains supported, but %s is the ' ...
     'recommended method for new work. Silence this with ' ...
     'warning(''off'',''RAVEN:legacyMethod'').'], oldFunction, char(preferred));
end
