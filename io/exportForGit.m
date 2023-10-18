function out=exportForGit(model,prefix,path,formats,mainBranchFlag,subDirs,cobraText)
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
%   mainBranchFlag      logical, if true, function will error if RAVEN (and
%                       COBRA if detected) is/are not on the main branch.
%                       (opt, default false)
%   subDirs             logical, whether model files for each file format 
%                       should be written in its own subdirectory, with
%                       'model' as parent directory, in accordance to the
%                       standard-GEM repository format. If false, all files
%                       are stored in the same folder. (opt, default true)
%   cobraText           logical, whether the txt file should be in COBRA
%                       Toolbox format using metabolite IDs, instead of
%                       metabolite names and compartments. (opt, default
%                       false)
%
%   Usage: exportForGit(model,prefix,path,formats,mainBranchFlag)
if nargin<7
    cobraText=false;
end
if nargin<6
    subDirs=true;
end
if nargin<5
    mainBranchFlag=false;
end
if nargin<4 || isempty(formats)
    formats={'mat', 'txt', 'xlsx', 'xml', 'yml'};
else
    formats=convertCharArray(formats);
end
if any(~ismember(formats, {'mat', 'txt', 'xlsx', 'xml', 'yml'}))
    EM='Unknown file format defined. Only mat, txt, xlsx, xml and yml are allowed file formats.';
    error(EM)
end
if nargin<3
    path='.';
else
    path=char(path);
end
if nargin<2
    prefix='model';
else
    prefix=char(prefix);
end

%Sort reactions, metabolites and genes alphabetically
model=sortIdentifiers(model);

%Get versions or commits of toolboxes:
RAVENver = getToolboxVersion('RAVEN','ravenCobraWrapper.m',mainBranchFlag);
COBRAver = getToolboxVersion('COBRA','initCobraToolbox.m',mainBranchFlag);

%Retrieve libSBML version:
[ravenDir,prevDir]=findRAVENroot();
cd(fullfile(ravenDir,'software','libSBML'));
try % 5.17.0 and newer
    libSBMLver=OutputSBML;
    libSBMLver=libSBMLver.libSBML_version_string;
catch % before 5.17.0
    fid = fopen('tempModelForLibSBMLversion.xml','w+');
    fclose(fid);
    evalc('[~,~,libSBMLver]=TranslateSBML(''tempModelForLibSBMLversion.xml'',0,0)');
    libSBMLver=libSBMLver.libSBML_version_string;
    delete('tempModelForLibSBMLversion.xml');
end
cd(prevDir)

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
    if cobraText==true
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
        exportModel(model,fullfile(filePath{5},strcat(prefix,'.xml')));
end

%Save file with versions:
fid = fopen(fullfile(path,'dependencies.txt'),'wt');
fprintf(fid,['MATLAB\t' version '\n']);
fprintf(fid,['libSBML\t' libSBMLver '\n']);
fprintf(fid,['RAVEN_toolbox\t' RAVENver '\n']);
if ~isempty(COBRAver)
    fprintf(fid,['COBRA_toolbox\t' COBRAver '\n']);
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
