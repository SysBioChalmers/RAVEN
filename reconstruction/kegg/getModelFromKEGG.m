function [model,KOModel]=getModelFromKEGG(varargin)
% getModelFromKEGG  Load the global KEGG model.
%
% Loads the global KEGG reaction/gene model from keggModel.mat. On first
% use --- when no keggModel.mat is present yet --- the underlying
% artefacts (the gene-free reference model plus the KO/reaction/organism-
% gene relational tables, published as a raven-data release) are
% downloaded and assembled instead, and the result is cached to
% keggModel.mat so later calls load instantly. The first build can take
% a while (the organism-gene table covers every KEGG organism) and needs
% a few hundred MB of disk.
%
% Name-Value Arguments
% --------------------
% keepSpontaneous : logical
%     include reactions labeled as "spontaneous" (default true).
% keepUndefinedStoich : logical
%     include reactions in the form n A <=> n+1 A (default true).
% keepIncomplete : logical
%     include reactions labelled as "incomplete", "erroneous" or "unclear"
%     (default true).
% keepGeneral : logical
%     include reactions labelled as "general reaction" (default false).
%
% Returns
% -------
% model : struct
%     a model structure with all reactions and the metabolites used in them.
% KOModel : struct
%     a model structure representing KEGG Orthology ids and their associated
%     genes. KO ids are stored as reactions.
%
% Notes
% -----
% The model output can be used as a template for fillGaps. In that case,
% remove the genes and rxnGeneMat fields before parsing:
%     model = rmfield(model, 'genes'), etc.
%
% See also
% --------
% getKEGGModelForOrganism, getPhylDist, fillGaps

ravenPath=findRAVENroot();

p=parseRAVENargs(varargin, {'keepSpontaneous',true; 'keepUndefinedStoich',true; ...
    'keepIncomplete',true; 'keepGeneral',false});
keepSpontaneous=p.keepSpontaneous;
keepUndefinedStoich=p.keepUndefinedStoich;
keepIncomplete=p.keepIncomplete;
keepGeneral=p.keepGeneral;

modelFile=fullfile(ravenPath,'reconstruction','kegg','keggModel.mat');
if isfile(modelFile)
    fprintf(['Importing the global KEGG model from ' strrep(modelFile,'\','/') '... ']);
    load(modelFile,'model','KOModel','isSpontaneous','isUndefinedStoich','isIncomplete','isGeneral');
    fprintf('COMPLETE\n');
else
    [model,KOModel,isSpontaneous,isUndefinedStoich,isIncomplete,isGeneral]=buildGlobalKEGGModel(ravenPath);
    fprintf(['Saving the global KEGG model to ' strrep(modelFile,'\','/') ' for future use... ']);
    save(modelFile,'model','KOModel','isSpontaneous','isUndefinedStoich','isIncomplete','isGeneral','-v7.3');
    fprintf('COMPLETE\n');
end

if keepSpontaneous==false
    model=removeReactions(model,intersect(isSpontaneous,model.rxns),true,true);
end
if keepUndefinedStoich==false
    model=removeReactions(model,intersect(isUndefinedStoich,model.rxns),true,true);
end
if keepIncomplete==false
    model=removeReactions(model,intersect(isIncomplete,model.rxns),true,true);
end
if keepGeneral==false
    model=removeReactions(model,intersect(isGeneral,model.rxns),true,true);
end
end
