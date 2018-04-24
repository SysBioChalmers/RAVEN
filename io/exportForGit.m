function out=exportForGit(model,prefix,path,formats)
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
%   formats             cell array of strings specifying in what file formats
%                       the model should be exported (opt, default to all
%                       formats as {'mat', 'txt', 'xlsx', 'xml', 'yml'})
%
%   Usage: exportForGit(model,prefix,path,formats)
%
%   Eduard Kerkhoven, 2018-04-24
%
if nargin<4
    formats={'mat', 'txt', 'xlsx', 'xml', 'yml'};
end
if ischar(formats)
    formats={formats};
end
if any(~ismember(formats, {'mat', 'txt', 'xlsx', 'xml', 'yml'}))
    EM='Unknown file format defined. Only mat, txt, xlsx, xml and yml are allowed file formats.';
    error(EM)
end
if nargin<3
    path='.';
end
if nargin<2
    prefix='model';
end

% Make folder structure if needed
folders={'ModelFiles','ComplementaryScripts','ComplementaryData'};
for i=1:length(folders);
    if ~exist(fullfile(path,folders{i}),'dir')
        mkdir(fullfile(path,folders{i}));
    end
end

for i=1:length(formats);
    if ~exist(fullfile(path,'ModelFiles',formats{i}),'dir')
        mkdir(fullfile(path,'ModelFiles',formats{i}))
    end
end

% Write MAT format
if ismember('mat', formats)
    save([fullfile(path,'ModelFiles','mat',prefix),'.mat'],'model');
end

% Write TXT format
if ismember('txt', formats)
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
end

% Write XLSX format
if ismember('xlsx', formats)
    exportToExcelFormat(model,strcat(prefix,'.xlsx'));
    movefile([prefix,'.xlsx'],fullfile(path,'ModelFiles','xlsx'));
end

% Write XML format
if ismember('xml', formats)
    exportModel(model,strcat(prefix,'.xml'));
    movefile([prefix,'.xml'],fullfile(path,'ModelFiles','xml'));
end

% Write YML format
if ismember('yml', formats)
    writeYaml(model,strcat(prefix,'.yml'));
    movefile([prefix,'.yml'],fullfile(path,'ModelFiles','yml'));
end

%Code below is modified from SysBioChalmers/yeast-GEM
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
%Retrieve libSBML version:
fid = fopen('tempModelForLibSBMLversion.xml','w+');
fclose(fid);
evalc('[~,~,libSBMLver]=TranslateSBML(''tempModelForLibSBMLversion.xml'',0,0)');
libSBMLver=libSBMLver.libSBML_version_string;
delete('tempModelForLibSBMLversion.xml');

%Save file with versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['MATLAB\t' version '\n']);
fprintf(fid,['libSBML\t' libSBMLver '\n']);
fprintf(fid,['RAVEN_toolbox\t' RAVENver '\n']);
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

movefile('*.txt',fullfile(path,'ModelFiles'));
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
