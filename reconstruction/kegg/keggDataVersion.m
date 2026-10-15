function version=keggDataVersion()
% keggDataVersion  The toolbox-assets release tag for KEGG artefacts.
%
% Reads reconstruction/kegg/KEGG_VERSION.md --- the single place this
% version is recorded --- so every function that downloads a KEGG
% artefact (getModelFromKEGG, getKEGGModelForOrganism) stays in sync
% when toolbox-assets publishes a new KEGG release. The version is the
% file's first line.
%
% Returns
% -------
% version : char
%     the toolbox-assets release tag, e.g. 'kegg118'.

ravenPath=findRAVENroot();
versionFile=fullfile(ravenPath,'reconstruction','kegg','KEGG_VERSION.md');
fid=fopen(versionFile,'r');
if fid==-1
    error('keggDataVersion:fileNotFound','Cannot read %s.',strrep(versionFile,'\','/'))
end
version=strtrim(fgetl(fid));
fclose(fid);
end
