function md5Hash=getMD5Hash(inputFile,varargin)
% getMD5Hash  Calculate the MD5 hash for a file.
%
% Parameters
% ----------
% inputFile : char
%     string with the path to the file for which the MD5 hash should be
%     calculated.
%
% Name-Value Arguments
% --------------------
% binEnd : char
%     string that indicates the operating system running on the client's
%     computer. Use ".exe" for Windows, ".mac" for macOS or leave it blank
%     for Linux (""). (default: the function automatically detects the
%     client's operating system).
%
% Returns
% -------
% md5Hash : char
%     string containing an MD5 hash for inputFile.
%
% Examples
% --------
%     md5Hash = getMD5Hash(inputFile, binEnd);
inputFile=char(inputFile);

p=parseRAVENargs(varargin, {'binEnd',[]});
binEnd=p.binEnd;
if isempty(binEnd)
    if isunix
        if ismac
            binEnd='.mac';
        else
            binEnd='';
        end
    elseif ispc
        binEnd='.exe';
    else
        error('Unknown OS, exiting.')
    end
else
    binEnd=char(binEnd);
end

%Check if binEnd is valid
if ~any(strcmp(binEnd,{'.mac','','.exe'}))
   error('Unknown OS, exiting.')
end

%Check file existence
inputFile=checkFileExistence(inputFile);

%Get string containing an MD5 hash
switch binEnd
    case '.mac'
        [~, md5HashMessage]=system(['md5 "' inputFile '"']);
    case ''
        [~, md5HashMessage]=system(['md5sum "' inputFile '"']);
    case '.exe'
        [~, md5HashMessage]=system(['certutil -hashfile "' inputFile '" MD5"']);
end

%Extract MD5 hash from a string
md5Hash = char(regexp(md5HashMessage,'[a-f0-9]{32}','match'));
end
