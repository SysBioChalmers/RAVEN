function localPath=fetchRavenDataAsset(destDir,releaseTag,assetName)
% fetchRavenDataAsset  Download (if needed) a raven-data release asset.
%
% Downloads assetName from the given raven-data release tag into destDir
% and returns its local path. An already-downloaded copy is reused as-is.
% Generic across raven-data's release families (KEGG artefacts, HMM
% libraries, binaries, ...): nothing here is specific to any one of
% them, so callers supply the release tag and asset name rather than
% this function hardcoding either.
%
% Parameters
% ----------
% destDir : char
%     directory to download into (created if it does not exist yet).
% releaseTag : char
%     the raven-data release tag, e.g. 'kegg118'.
% assetName : char
%     the release asset file name, e.g. 'kegg118_core.tar.gz'.
%
% Returns
% -------
% localPath : char
%     path to the downloaded file.

if ~isfolder(destDir)
    mkdir(destDir);
end
localPath=fullfile(destDir,assetName);
if isfile(localPath)
    return;
end
fprintf(['Downloading ' assetName '... ']);
try
    websave(localPath,['https://github.com/SysBioChalmers/raven-data/releases/download/' releaseTag '/' assetName]);
catch ME
    if strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError')
        error('Failed to download %s, the server returned a 404 error, try again later. If the problem persists please report it on the RAVEN GitHub Issues page: https://github.com/SysBioChalmers/RAVEN/issues',assetName)
    else
        rethrow(ME);
    end
end
fprintf('COMPLETE\n');
end
