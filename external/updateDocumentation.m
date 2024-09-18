function updateDocumentation()
% updateDocumentation
%	Updates HTML documentation files for all RAVEN functions. It should
%	only be used when preparing `develop` branch for a new RAVEN release
%
%	Usage: updateDocumentation()

%Get the RAVEN path
ravenDir=findRAVENroot();
%Make sure that RAVEN-provided m2html is used
path(fullfile(ravenDir,'software','m2html'),path);
%Get a non-redundant list of RAVEN subdirectories containing MATLAB
%functions. Absolute paths are not compatible with M2HTML, so convert them
%to the relative paths instead.
ravenDirs=dir(fullfile(ravenDir,'**/*.m'));
ravenDirs=unique({ravenDirs.folder})';

%Get rid of MATLAB functions from external software
ravenDirs(startsWith(ravenDirs,strcat(ravenDir,filesep,'software')))=[];
ravenDirs(startsWith(ravenDirs,strcat(ravenDir,filesep,'legacy',filesep,'software')))=[];

%Remove keggModel.mat if it exists
if exist(fullfile(ravenDir,'external','kegg','keggModel.mat'), 'file') == 2
    delete(fullfile(ravenDir,'external','kegg','keggModel.mat'));
end

%Remove existing "doc" directory from RAVEN
rmdir(fullfile(ravenDir,'doc'),'s');

%Make relative path
relStart = numel(ravenDir)+2;
for i=1:numel(ravenDirs)
    ravenDirs{i,1} = ravenDirs{i,1}(relStart:end);
end

%Save the current working directory and go to RAVEN root directory
originalDir=pwd;
cd(ravenDir);
%Generate HTML documentation files for RAVEN MATLAB functions
m2html('mFiles',ravenDirs,'htmldir','doc');
%Go back to the original working directory
cd(originalDir);

end
