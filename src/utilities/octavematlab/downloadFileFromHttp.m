function downloadFileFromHttp(url,fileLocation)
% downloadFileFromHttp
%   Dowload a file from a http url and store at the indicated location.
%
%   Usage: downloadFileFromHttp(url,fileLocation)

if isOctave
    urlwrite(url,fileLocation);
else
    websave(fileLocation,url);
end
end
    