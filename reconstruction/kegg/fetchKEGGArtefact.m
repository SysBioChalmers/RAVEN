function localPath=fetchKEGGArtefact(ravenPath,assetName)
% fetchKEGGArtefact  Download (if needed) a raven-data kegg118 release asset.
%
% Downloads assetName from the kegg118 raven-data release into
% <ravenPath>/reconstruction/kegg/ and returns its local path. An
% already-downloaded copy is reused as-is.
%
% Parameters
% ----------
% ravenPath : char
%     the RAVEN root directory, as returned by findRAVENroot.
% assetName : char
%     the raven-data release asset file name, e.g. 'kegg118_core.tar.gz'.
%
% Returns
% -------
% localPath : char
%     path to the downloaded file.

localPath=fullfile(ravenPath,'reconstruction','kegg',assetName);
if isfile(localPath)
    return;
end
fprintf(['Downloading ' assetName '... ']);
try
    websave(localPath,['https://github.com/SysBioChalmers/raven-data/releases/download/kegg118/' assetName]);
catch ME
    if strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError')
        error('Failed to download %s, the server returned a 404 error, try again later. If the problem persists please report it on the RAVEN GitHub Issues page: https://github.com/SysBioChalmers/RAVEN/issues',assetName)
    else
        rethrow(ME);
    end
end
fprintf('COMPLETE\n');
end
