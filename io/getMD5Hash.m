function md5Hash=getMD5Hash(inputFile,binEnd)
% getMD5Hash
%   Calculates MD5 hash for a file
%
%   Input:
%   inputFile       string with the path to file for which MD5 hash should
%                   be calculated
%   binEnd          string that shows the operating system running in the
%                   client's computer. Use ".exe" for Windows, ".mac" for
%                   macOS or leave it blank for Linux (""). (opt, by
%                   default the function automatically detects the client's
%                   operating system)
%
%   Output:
%   md5Hash         string containing an MD5 hash for inputFile
%   
%   Usage: md5Hash=getMD5Hash(inputFile,binEnd)
inputFile=char(inputFile);

if nargin<2
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
