function updateDocumentation()
% updateDocumentation
%	Regenerates the HTML documentation in RAVEN's doc directory from the
%	function help texts. The doc directory is deleted and rebuilt, so
%	anything hand-edited there is lost.
%
%	The Documentation workflow runs this on every push to a release branch
%	and commits the result, so doc stays in step with the sources without
%	anyone running it by hand. Run it locally only to preview.
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

%Skip deprecated wrappers: they are on the path so existing scripts keep
%working, but documenting them would advertise functions that are on their
%way out. Each one names its replacement when called.
ravenDirs(startsWith(ravenDirs,strcat(ravenDir,filesep,'deprecated')))=[];

%Remove keggModel.mat if it exists
if exist(fullfile(ravenDir,'reconstruction','kegg','keggModel.mat'), 'file') == 2
    delete(fullfile(ravenDir,'reconstruction','kegg','keggModel.mat'));
end

%Remove the existing "doc" directory from RAVEN
docDir=fullfile(ravenDir,'doc');
if isfolder(docDir)
    rmdir(docDir,'s');
end

%Make relative path
relStart = numel(ravenDir)+2;
for i=1:numel(ravenDirs)
    ravenDirs{i,1} = ravenDirs{i,1}(relStart:end);
end

%Recreate the output tree before handing over to m2html. m2html only creates an
%output directory when exist() says it is not already there, and exist() keeps
%answering "directory" for a relative path that was just deleted -- a stale
%answer that neither rmpath nor rehash clears. It would therefore skip the
%mkdir and then fail writing the first file into a directory that is gone.
mkdir(docDir);
for i=1:numel(ravenDirs)
    if ~isempty(ravenDirs{i})
        mkdir(fullfile(docDir,ravenDirs{i}));
    end
end

%Save the current working directory and go to RAVEN root directory
originalDir=pwd;
cd(ravenDir);
%Generate HTML documentation files for RAVEN MATLAB functions
m2html('mFiles',ravenDirs,'htmldir','doc');
%Go back to the original working directory
cd(originalDir);

end
