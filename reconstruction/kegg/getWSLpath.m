function path=getWSLpath(path)
% getWSLpath  Translate a Windows-style path to its Unix WSL equivalent.
%
% Translate a Windows-style path to its Unix WSL (Windows Subsystem for
% Linux) equivalent.
%
% Parameters
% ----------
% path : char
%     string with directory of file path, in Windows-style (e.g.
%     'C:\Directory\').
%
% Returns
% -------
% path : char
%     string with directory of file path, in Unix style (e.g.
%     '/mnt/c/Directory/').
%
% Examples
% --------
%     path = getWSLpath(path);
%
% Notes
% -----
% Uses the WSL function 'wslpath' to translate the path.
[status,path]=system(['wsl wslpath ''' path '''']);
if status==-1
    error('Cannot get access to Windows Subsystem for Linux, check your WSL installation')
end
path=path(1:end-1);% Remove final character (line-break)
end
