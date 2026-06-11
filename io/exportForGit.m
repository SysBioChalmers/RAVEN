function out=exportForGit(model,varargin)
% exportForGit  Export a model for a Git-maintained model repository.
%
% Generates a directory structure and populates it with model files, ready
% to be committed to a Git(Hub) maintained model repository. Writes the
% model as SBML L3V1 FBCv2 (both XML and YAML), COBRA text, Matlab MAT-file
% and Microsoft Excel formats.
%
% Parameters
% ----------
% model : struct
%     model structure in RAVEN format that should be exported.
%
% Name-Value Arguments
% --------------------
% prefix : char
%     prefix for all filenames (default 'model').
% path : char
%     path where the directory structure should be generated and populated
%     with all files (default current working directory).
% formats : cell
%     cell array of strings specifying in what file formats the model
%     should be exported (default all formats as {'mat', 'txt', 'xlsx',
%     'xml', 'yml'}).
% mainBranchFlag : logical
%     if true, function will error if RAVEN (and GECKO for enzyme-
%     constrained models) is/are not on the main branch (default false).
% subDirs : logical
%     whether model files for each file format should be written in their
%     own subdirectory, with 'model' as parent directory, in accordance to
%     the standard-GEM repository format. If false, all files are stored in
%     the same folder (default true).
% COBRAtext : logical
%     whether the txt file should be in COBRA Toolbox format using
%     metabolite IDs, instead of metabolite names and compartments
%     (default false).
% neverPrefixIDs : logical
%     true if prefixes are never added to identifiers, even if they start
%     with e.g. digits. This might result in invalid SBML files (default
%     false).
%
% Examples
% --------
%     exportForGit(model, prefix, path, formats, mainBranchFlag, subDirs, ...
%         COBRAtext, neverPrefixIDs);
p=parseRAVENargs(varargin, {'prefix',[]; 'path',[]; 'formats',[]; 'mainBranchFlag',[]; 'subDirs',[]; 'COBRAtext',[]; 'neverPrefixIDs',false});
prefix=p.prefix; path=p.path; formats=p.formats; mainBranchFlag=p.mainBranchFlag; subDirs=p.subDirs; COBRAtext=p.COBRAtext; neverPrefixIDs=p.neverPrefixIDs;
if isempty(COBRAtext)
    COBRAtext=false;
end
if isempty(subDirs)
    subDirs=true;
end
if isempty(mainBranchFlag)
    mainBranchFlag=false;
end
if isempty(formats)
    formats={'mat', 'txt', 'xlsx', 'xml', 'yml'};
else
    formats=convertCharArray(formats);
end
if any(~ismember(formats, {'mat', 'txt', 'xlsx', 'xml', 'yml'}))
    EM='Unknown file format defined. Only mat, txt, xlsx, xml and yml are allowed file formats.';
    error(EM)
end
if isempty(path)
    path='.';
else
    path=char(path);
end
if isempty(prefix)
    prefix='model';
else
    prefix=char(prefix);
end

%Sort reactions, metabolites and genes alphabetically
model=sortIdentifiers(model);

%Get versions or commits of toolboxes. GECKO is only relevant for
%enzyme-constrained models (those with an "ec" field).
RAVENver = getToolboxVersion('RAVEN','ravenCobraWrapper.m',mainBranchFlag);
if isfield(model,'ec')
    GECKOver = getToolboxVersion('GECKO','GECKOInstaller.m',mainBranchFlag);
else
    GECKOver = [];
end
%Retrieve libSBML version:
[ravenDir,prevDir]=findRAVENroot();
try % 5.17.0 and newer
    libSBMLver=OutputSBML_RAVEN;
    libSBMLver=libSBMLver.libSBML_version_string;
catch % before 5.17.0
    fid = fopen('tempModelForLibSBMLversion.xml','w+');
    fclose(fid);
    evalc('[~,~,libSBMLver]=TranslateSBML_RAVEN(''tempModelForLibSBMLversion.xml'',0,0)');
    libSBMLver=libSBMLver.libSBML_version_string;
    delete('tempModelForLibSBMLversion.xml');
end

% Make models folder, no warnings if folder already exists
if subDirs
    path=fullfile(path,'model');
    filePath=strcat(path,filesep,{'txt','yml','mat','xlsx','xml'});
    [~,~,~]=mkdir(path);
    for i = 1:length(formats)
        [~,~,~]=mkdir(fullfile(path,formats{i}));
    end
else
    filePath=cell(1,5); filePath(:)={path};
end


% Write TXT format
if ismember('txt', formats)
    fid=fopen(fullfile(filePath{1},strcat(prefix,'.txt')),'w');
    if COBRAtext==true
        eqns=constructEquations(model,model.rxns,false,false,false);
        eqns=strrep(eqns,' => ','  -> ');
        eqns=strrep(eqns,' <=> ','  <=> ');
        eqns=regexprep(eqns,'> $','>');
        grRules=regexprep(model.grRules,'\((?!\()','( ');
        grRules=regexprep(grRules,'(?<!\))\)',' )');
    else
        eqns=constructEquations(model,model.rxns);
        grRules=model.grRules;
    end
    fprintf(fid, 'Rxn name\tFormula\tGene-reaction association\tLB\tUB\tObjective\n');
    for i = 1:numel(model.rxns)
        fprintf(fid, '%s\t', model.rxns{i});
        fprintf(fid, '%s \t', eqns{i});
        fprintf(fid, '%s\t', grRules{i});
        fprintf(fid, '%6.2f\t%6.2f\t%6.2f\n', model.lb(i), model.ub(i), model.c(i));
    end
    fclose(fid);
end

% Write YML format
if ismember('yml', formats)
    writeYAMLmodel(model,fullfile(filePath{2},strcat(prefix,'.yml')));
end

% Write MAT format
if ismember('mat', formats)
    save(fullfile(filePath{3},strcat(prefix,'.mat')),'model');
end

% Write XLSX format
if ismember('xlsx', formats)
    exportToExcelFormat(model,fullfile(filePath{4},strcat(prefix,'.xlsx')));
end

% Write XML format
if ismember('xml', formats)
        exportModel(model,fullfile(filePath{5},strcat(prefix,'.xml')),neverPrefixIDs);
end

%Save file with versions:
fid = fopen(fullfile(path,'dependencies.txt'),'wt');
fprintf(fid,['MATLAB\t' version '\n']);
fprintf(fid,['libSBML\t' libSBMLver '\n']);
fprintf(fid,['RAVEN_toolbox\t' RAVENver '\n']);
if ~isempty(GECKOver)
    fprintf(fid,['GECKO_toolbox\t' GECKOver '\n']);
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

function version = getToolboxVersion(toolbox,fileID,mainBranchFlag)
% getToolboxVersion  Return the version of RAVEN or GECKO.
%
% Returns the version of the toolbox (from its version.txt), or if not
% available the latest commit hash (7 characters).
%
% toolbox        name of the toolbox (e.g. 'RAVEN' or 'GECKO').
% fileID         a file only found in that toolbox, used to locate its root
%                (e.g. 'ravenCobraWrapper.m' or 'GECKOInstaller.m').
% mainBranchFlag if true, error if the toolbox is not on the main branch.
currentPath = pwd;
version     = '';

%Try to find root of toolbox:
try
    toolboxPath = which(fileID);                %full file path
    slashPos    = getSlashPos(toolboxPath);
    toolboxPath = toolboxPath(1:slashPos(end)); %folder path
    %Go up until the root is found:
    D = dir(toolboxPath);
    while ~ismember({'.git'},{D.name})
        slashPos    = getSlashPos(toolboxPath);
        toolboxPath = toolboxPath(1:slashPos(end-1));
        D           = dir(toolboxPath);
    end
    cd(toolboxPath);
catch
    disp([toolbox ' toolbox cannot be found'])
    version = 'unknown';
end
%Check if in main:
if mainBranchFlag
    [~,currentBranch] = system('git rev-parse --abbrev-ref HEAD');
    currentBranch = strtrim(currentBranch);
    if any([strcmp(currentBranch, "main"), strcmp(currentBranch, "master")])
        cd(currentPath);
        error(['ERROR: ' toolbox ' not in main (or master) branch. Check-out this branch of ' toolbox ' before submitting model for Git.'])
    end
end
%Try to find version file of the toolbox:
if isempty(version)
    try
        fid     = fopen([toolboxPath 'version.txt'],'r');
        version = fscanf(fid,'%s');
        fclose(fid);
    catch
        %If no file available, look up the tag:
        try
            [~,version] = system('git describe --tags');
            version = strtrim(version);
            [~,commit] = system('git log -n 1 --format=%H');
            commit = commit(1:7);
            %If no tag available or commit is part of tag, get commit instead:
            if ~isempty(strfind(version,'fatal')) || ~isempty(strfind(version,commit))
                version = ['commit ' commit];
            else
                version = strrep(version,'v','');
            end
        catch
            version = 'unknown';
        end
    end
end
cd(currentPath);
end

function slashPos = getSlashPos(path)
slashPos = strfind(path,'\');       %Windows
if isempty(slashPos)
    slashPos = strfind(path,'/');   %MAC/Linux
end
end
