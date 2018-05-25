function addRavenToUserPath(overwrite)
% This function writes a startup.m file in the userpath, adding RAVEN (and
% all subdirectories) to the path each time Matlab is started.
% This function is useful if the user has no rights to save paths to the
% pathdef.m file in the matlabroot, for instance due to a lack of admin
% rights. As the startup.m file in the userpath automatically runs with
% each Matlab start, the paths are automatically loaded.
%
%   overwrite       logical, whether startup.m in the userpath should
%                   overwritten (otherwise the RAVEN paths are appended)
%                   (opt, default true)
%
%   Eduard Kerkhoven, 2018-04-26
%

if nargin1
    overwrite=true;
end

% Get current RAVEN directory
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

% Lists all subdirectories
subpath=regexp(genpath(ravenDir),pathsep,'split');
% Remove .git and doc folders
pathsToKeep=cellfun(@(x) isempty(strfind(x,'.git')),subpath) & cellfun(@(x) isempty(strfind(x,'doc')),subpath);
% Only keep useful paths
subpath = subpath(pathsToKeep);
subpath = subpath(1:end-1); % Remove last entry, is empty field

% Write startup.m file
if overwrite
    fid=fopen(fullfile(userpath,'startup.m'),'w');
    fprintf(fid,'%sn','%%% RAVEN path');
else
    fid=fopen(fullfile(userpath,'startup.m'),'a');
    fprintf(fid,'n%sn','%%% RAVEN path');
end
fprintf(fid,'%sn',strcat('addpath(''',subpath{1},''',...'));
for i=2(length(subpath)-1)
    fprintf(fid,'t%sn',strcat('''',subpath{i},''',...'));
end
fprintf(fid,'t%s',strcat('''',subpath{length(subpath)},''');'));
fclose(fid);
end
