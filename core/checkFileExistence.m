function files=checkFileExistence(files,makeFullPath,allowSpace)
% checkFileExistence
%   Check whether files exist. If no full path is given a file should be
%   located in the current folder, which by default is appended to the
%   filename.
%
%   Input:
%   files           string or cell array of strings with path to file(s) or
%                   path or filename(s)
%   makeFullPath    logical, whether files located in the current folder
%                   should be provided with the full path (opt, default
%                   true)
%   allowSpace      logical, whether 'space' character is allowed in the
%                   path (opt, default true)
%   
%   Usage: files=checkFileExistence(files,makeFullPath,allowSpace)
%
%   Eduard Kerkhoven, 2019-09-30
%

if nargin<3
    allowSpace = true;
end
if nargin<2
    makeFullPath = true;
end
if isstr(files)
    files={files};
end
filesOriginal = files;

inCurrDir = ~contains(files,'\') & ~contains(files,'/');
files(inCurrDir) = fullfile(cd,files(inCurrDir));

for i=1:numel(files)
    if ~exist(files{i},'file')
        error('File "%s" cannot be found\n',files{i});
    elseif allowSpace == true & strfind(files{i},' ')
        error('File "%s" has a space in the filename or path. Remove this before running getDiamond\n',files{i});
    end
end

if makeFullPath == false;
    files = filesOriginal;
end