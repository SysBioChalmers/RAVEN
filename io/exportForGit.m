function out=exportForGit(model,prefix,path)

if nargin<3
    path='.';
end
if nargin<2
    prefix='model';
end

%% Make folder structure if needed
if ~exist(fullfile(path,'ModelFiles'),'dir')
    mkdir(fullfile(path,'ModelFiles'));
end

formats={'xml','yaml','txt','mat'};

for i=1:length(formats);
    if ~exist(fullfile(path,'ModelFiles',formats{i}),'dir')
        mkdir(fullfile(path,'ModelFiles',formats{i}))
    end
end

%% Write txt format
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

%% Write XML (SBML) and YAML formats
exportModel(model,prefix,true);
movefile([prefix,'.xml'],fullfile(path,'ModelFiles','xml'));
movefile([prefix,'.yml'],fullfile(path,'ModelFiles','yaml'));
save([fullfile(path,'ModelFiles','mat',prefix),'.mat'],'model');