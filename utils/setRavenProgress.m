function setRavenProgress(mode)
% setRavenProgress  Set the RAVEN progress-reporting backend.
%
% Sets the preference that progressReport uses to decide how progress is
% shown. The setting is persistent across MATLAB sessions.
%
% Parameters
% ----------
% mode : char
%     one of:
%
%     - 'auto' : choose automatically (default). GUI on an interactive
%       desktop, an animated bar in an interactive terminal, plain log lines
%       on a headless / -batch server, and no output during unit tests.
%     - 'gui'  : always use a MATLAB waitbar.
%     - 'cli'  : always use the animated command-line bar.
%     - 'log'  : always print plain milestone lines (safe in redirected
%       logs and on headless servers).
%     - 'none' / 'off' : never show progress.
%
% Examples
% --------
%     setRavenProgress('log');   % e.g. on a compute cluster
%     setRavenProgress('off');   % silence all progress reporting
%     setRavenProgress('auto');  % restore the default behaviour
%
% See also
% --------
% progressReport

if nargin < 1 || isempty(mode)
    mode = 'auto';
end
mode = lower(char(mode));
if strcmp(mode,'off')
    mode = 'none';
end
if ~ismember(mode,{'auto','gui','cli','log','none'})
    error('Unknown progress mode "%s". Use auto, gui, cli, log or off.', mode);
end
setpref('RAVEN','progressBar',mode);
end
