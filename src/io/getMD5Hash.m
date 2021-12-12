function md5Hash=getMD5Hash(inputFile)
% getMD5Hash
%   Calculates MD5 hash for a file
%
%   Input:
%   inputFile       string with the path to file for which MD5 hash should
%                   be calculated
%
%   Output:
%   md5Hash         string containing an MD5 hash for inputFile
%   
%   Usage: md5Hash=getMD5Hash(inputFile)

%Get the OS specific binary ending (e.g. exe for Windows)
binEnd = binaryEnding();

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
