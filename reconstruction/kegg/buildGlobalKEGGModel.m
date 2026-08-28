function [model,KOModel,isSpontaneous,isUndefinedStoich,isIncomplete,isGeneral]=buildGlobalKEGGModel(ravenPath)
% buildGlobalKEGGModel  Download and assemble the global KEGG model.
%
% Fetches the raven-data <kegg version>_core.tar.gz bundle (the
% gene-free reference model plus the ko_reaction / organism_gene_ko /
% rxn_flags relational tables --- see raven-toolbox's
% docs/maintenance/kegg_data_format.md), then joins every organism's
% genes onto the reference reactions through their shared KO ids. This
% is the data getModelFromKEGG used to load from a pre-built
% keggModel.mat; see RAVEN issue #704. The KEGG version is read from
% keggDataVersion, not hardcoded here.
%
% Parameters
% ----------
% ravenPath : char
%     the RAVEN root directory, as returned by findRAVENroot.
%
% Returns
% -------
% model : struct
%     the full global KEGG model (all reactions/metabolites, genes and
%     rxnGeneMat spanning every KEGG organism). Callers narrow this down
%     to one organism (see getKEGGModelForOrganism).
% KOModel : struct
%     a minimal model struct whose rxns are the KO ids linked to at
%     least one reaction (getKEGGModelForOrganism's HMM-search path uses
%     this as a KO id lookup table; nothing else of KOModel is read).
% isSpontaneous, isUndefinedStoich, isIncomplete, isGeneral : cell arrays
%     reaction ids carrying the corresponding rxn_flags quality flag.

kver=keggDataVersion();
keggDir=fullfile(ravenPath,'reconstruction','kegg');
archive=fetchRavenDataAsset(keggDir,kver,[kver '_core.tar.gz']);

coreDir=fullfile(keggDir,[kver '_core']);
if ~isfolder(coreDir)
    fprintf('Extracting the KEGG core artefacts... ');
    untar(archive,coreDir);
    fprintf('COMPLETE\n');
end

fprintf('Reading the KEGG reference model... ');
refFileGz=fullfile(coreDir,[kver '_reference_model.yml.gz']);
refFile=refFileGz(1:end-3);
if ~isfile(refFile)
    gunzip(refFileGz);
end
model=readYAMLmodel(refFile);
fprintf('COMPLETE\n');

fprintf('Reading the KEGG relational tables... ');
koReaction=readKEGGTable(fullfile(coreDir,[kver '_ko_reaction.tsv.gz']));
rxnFlags=readKEGGTable(fullfile(coreDir,[kver '_rxn_flags.tsv.gz']));
organismGeneKO=readKEGGTable(fullfile(coreDir,[kver '_organism_gene_ko.tsv.gz']));
fprintf('COMPLETE\n');

isSpontaneous=flaggedReactions(rxnFlags,'spontaneous');
isUndefinedStoich=flaggedReactions(rxnFlags,'undefined_stoich');
isIncomplete=flaggedReactions(rxnFlags,'incomplete');
isGeneral=flaggedReactions(rxnFlags,'general');

KOModel.id='KOModel';
KOModel.description='KEGG Orthology ids linked to at least one reaction';
KOModel.rxns=unique(koReaction.ko);

fprintf('Joining organism genes onto the reference reactions (this can take a while for the full KEGG gene set)... ');
[model.genes,model.rxnGeneMat]=buildGlobalGPR(model.rxns,koReaction,organismGeneKO);
fprintf('COMPLETE\n');
end

function ids=flaggedReactions(rxnFlags,column)
mask=strcmpi(rxnFlags.(column),'true') | strcmp(rxnFlags.(column),'1');
ids=rxnFlags.reaction(mask);
end
