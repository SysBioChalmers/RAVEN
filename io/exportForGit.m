function out=exportForGit(model,prefix,path)
% exportForGit
%   Generates a directory structure and populates this with model files, ready
%   to be commited to a Git(Hub) maintained model repository. Writes the model
%   as SBML L3V1 FBCv2 (both XML and YAML), COBRA text, Matlab MAT-file
%   orthologies in KEGG
%
%   model               model structure in RAVEN format that should be exported
%   prefix              prefix for all filenames (opt, default 'model')
%   path                path where the directory structure should be generated
%                       and populated with all files (opt, default to current
%                       working directory)
%
%   Usage: exportForGit(model,prefix,path)
%
%   Eduard Kerkhoven, 2018-03-14
%
if nargin<3
    path='.';
end
if nargin<2
    prefix='model';
end

% Make folder structure if needed
if ~exist(fullfile(path,'ModelFiles'),'dir')
    mkdir(fullfile(path,'ModelFiles'));
end

formats={'xml','yaml','txt','mat'};

for i=1:length(formats);
    if ~exist(fullfile(path,'ModelFiles',formats{i}),'dir')
        mkdir(fullfile(path,'ModelFiles',formats{i}))
    end
end

% Write txt format
fid=fopen([fullfile(path,'ModelFiles','txt',prefix),'.txt'],'w');
eqns=constructEquations(model,model.rxns,false,false,false,true);
eqns=strrep(eqns,' => ','  -> ');
eqns=strrep(eqns,' <=> ','  <=> ');
eqns=regexprep(eqns,'> $','>');
grRules=regexprep(model.grRules,'\((?!\()','( ');
grRules=regexprep(grRules,'(?<!\))\)',' )');
fprintf(fid, 'Rxn name\tFormula\tGene-reaction association\tLB\tUB\tObjective\n');
for i = 1:numel(model.rxns)
    fprintf(fid, '%s\t', model.rxns{i});
    fprintf(fid, '%s \t', eqns{i});
    fprintf(fid, '%s\t', grRules{i});
    fprintf(fid, '%6.2f\t%6.2f\t%6.2f\n', model.lb(i), model.ub(i), model.c(i));
end
fclose(fid);

% Write XML (SBML) and YAML formats
exportModel(model,prefix,true);
movefile([prefix,'.xml'],fullfile(path,'ModelFiles','xml'));
movefile([prefix,'.yml'],fullfile(path,'ModelFiles','yaml'));
save([fullfile(path,'ModelFiles','mat',prefix),'.mat'],'model');

%Code below is modified from SysBioChalmers/YeastMetabolicNetwork-GEM
%Detect boundary metabolites and save them in a .txt file:
fid = fopen('boundaryMets.txt','wt');
for i = 1:length(model.rxns)
    pos = find(model.S(:,i) ~= 0);
    if length(pos) == 1 %Exchange rxn
        fprintf(fid,[model.mets{pos} '\t' model.metNames{pos} '\n']);
    end
end
fclose(fid);

%Track versions
RAVENver = getVersion('checkInstallation.m','version.txt');
%Retrieve latest COBRA commit:
COBRApath   = which('initCobraToolbox.m');
if ~isempty(COBRApath)
    slashPos    = getSlashPos(COBRApath);
    COBRApath   = COBRApath(1:slashPos(end)-1);
    currentPath = pwd;
    cd(COBRApath)
    try
        COBRAcommit = git('log -n 1 --format=%H');
    catch
        disp('COBRA is not fully installed (including Git wrapper)')
        COBRAcommit = 'unknown';
    end
    cd(currentPath)
else
    disp('COBRA version cannot be found')
end

%Save file with versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['RAVEN_toolbox\tv' RAVENver '\n']);
if ~isempty(COBRApath)
    fprintf(fid,['COBRA_toolbox\tcommit ' COBRAcommit(1:7) '\n']);
end
if isfield(model,'modelVersion')
    fields = fieldnames(model.modelVersion);
    for i = 1:length(fields)
        value = model.modelVersion.(fields{i});
        fprintf(fid,[fields{i} '\t' num2str(value) '\n']);
    end
end
fclose(fid);
end

function version = getVersion(IDfileName,VERfileName)
    try
        path     = which(IDfileName);
        slashPos = getSlashPos(path);
        path     = path(1:slashPos(end-1));
        fid      = fopen([path VERfileName],'r');
        version  = fscanf(fid,'%s');
        fclose(fid);
        catch
        version = '?';
    end
end

function slashPos = getSlashPos(path)
    slashPos = strfind(path,'\');       %Windows
    if isempty(slashPos)
        slashPos = strfind(path,'/');   %MAC/Linux
    end
end