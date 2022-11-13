% tutorial6
%   This exercise demonstrates how to reconstruct a combined draft GEM by
%   from KEGG and MetaCyc pathway databases. A combined model with the
%   comprehensive coverage of metabolic pathways is generated from
%   different de novo reconstruction approaches. The input is a FASTA
%   format file with whole-proteome sequences. The combined model is
%   subsequently used for refinement of existing high-quality model and
%   generation of a new version of GEM, by utilizing the manual curation
%   results. This tutorial is a showcase of the new features released in
%   RAVEN 2.0 through demonstrating the utilization of the newly developed
%   functions on GEM reconstruction and curation for Streptomyces
%   coelicolor strain A3(2). Users may apply this script as the template in
%   their own work for other organisms.
%   This refers to Tutorial 6 from "RAVEN tutorials.docx"

%Before reconstruction, a FASTA file with protein sequences of the target
%organism needs to be prepared. In this tutorial, all protein sequences of
%S. coelicolor A3(2) were downloaded from the NCBI genome database and
%provided in current folder (Sco_all_protein.faa). The description lines
%(starting with ">" character) of this FASTA file were modified by keeping
%only the locus tag (e.g. SCO6005), which is commonly used as the gene
%names in GEMs.

%Generate a draft GEM from MetaCyc database. The first parameter is
%organismID and needs to be specified by user. The other parameters are set
%to exclude unbalanced and undetermined reactions, but keep transport
%reactions. The two parameters for homology search are set to the default
%values that have been optimized to capture the protein hits with both
%comprehensive coverage and the least false positives.
ScoMetaCycDraftModel=getMetaCycModelForOrganism('ScoMetaCyc','Sco_all_protein.faa',1);

%Generate the two draft models from KEGG database. The first one is
%reconstructed from the genome information annotated by KEGG. Since KEGG
%provides genome annotations for over 5000 organisms, their GEMs can thus
%be generated without providing protein sequences while only using the
%organismID. The organisms with KEGG annotation and their associated ids
%(e.g. 'sco' for Streptomyces coelicolor) can be found from here:
%https://www.kegg.jp/kegg/catalog/org_list.html.
%Generate a draft model using KEGG annotation and exclude incomplete
%reactions and reactions with undefined stoichiometry
ScoKEGGAnnotation=getKEGGModelForOrganism('sco','','','',1,0,0);

%Generate the second KEGG-based draft model based on sequence homology to
%KEGG Ortholog sequence clusters while excluding incomplete reactions and
%the ones with undefined stoichiometry. Type "help getKEGGModelForOrganism"
%to see the detailed instructions for the choice of different parameters.
%The default values for homology search are used because they have been
%optimized for the best performance.
ScoKEGGHomology=getKEGGModelForOrganism('ScoKEGGHMMs','Sco_all_protein.faa','prok90_kegg102','',1,0,0);

%De novo reconstruction from MetaCyc should take about 10 minutes, while
%both reconstructions from KEGG may take up to 50-60 minutes in overall
%depending on the computer hardware and the system used. Given that KEGG
%and MetaCyc databases are formulated in different ways; KEGG relies on
%sequence-based annotation, while MetaCyc collects only experimentally
%verified pathways. Therefore, integration of MetaCyc- and KEGG-based draft
%models could have a better coverage of the metbolism for the target
%organism.

%At first, the two KEGG-based models can be directly merged
ScoKEGGDraftModel=mergeModels({ScoKEGGAnnotation ScoKEGGHomology});

disp(numel(ScoKEGGHomology.rxns)+numel(ScoKEGGAnnotation.rxns));
disp(numel(ScoKEGGDraftModel.rxns));
%By checking the reaction number, it can be seen that reaction number in
%the merged model equals adding up the reaction numbers in homology and
%annotation KEGG draft models. And there are duplicated reactions in this
%merged model.

%Merge the duplicated reactions into one and combine multiple iso-enzymes
%into a single grRules. The expandModel and contractModel functions are
%utilised, see their instructions for details.
ScoKEGGDraftModel=expandModel(ScoKEGGDraftModel);
ScoKEGGDraftModel=contractModel(ScoKEGGDraftModel);

%In the end, KEGG- and MetaCyc-based draft models can be further combined
%into an integrated GEM. This step is achieved by the function
%combineMetaCycKEGGModels, which firstly converts metabolite and reaction
%identifiers in the KEGG model into corresponding MetaCyc ids, and then
%detects duplications and keeps only unique reactions and metabolites that
%are mostly in MetaCyc identifiers.
ScoCombinedDraftModel=combineMetaCycKEGGModels(ScoMetaCycDraftModel, ScoKEGGDraftModel);

%With this combined model, the existing iMK1208 model is refined by
%incorporating the new pathways/reactions that are found in the combined
%model but absent from the previous iMK1208 model. At the time of
%publication, a total of 398 reactions in the combined draft model were
%determined as new pathways based on manual curation results, which have
%been organized into TableS3. Load these manually selected reactions and
%their subSystems as the pre-processed array structure selectedNewRxns.
load('iMK1208+suppInfo.mat','selectedNewRxns');
disp(selectedNewRxns);

%Generate a sub-model that includes only these new reactions. This is
%implemented by subtracting the other reactions from the combined model
%using function removeReactions
rxnsToRemove=setdiff(ScoCombinedDraftModel.rxns,selectedNewRxns.rxns);
newRxnSubModel=removeReactions(ScoCombinedDraftModel,rxnsToRemove,1,1);

%Now this newRxnSubModel contains only these new reactions. It needs to be
%modified before merging with the iMK1208. Since these reactions are
%metabolic ones and can all be assigned to the cytoplasm compartment
newRxnSubModel.comps{1}='c';
newRxnSubModel.compNames{1}='Cytoplasm';

%Amend to the genes and rxnGeneMat fields
newRxnSubModel.genes={};
for k=1:numel(newRxnSubModel.rxns)
   newRxnSubModel.genes=[newRxnSubModel.genes;transpose(strsplit(newRxnSubModel.grRules{k},' or '))];
end
newRxnSubModel.genes=unique(newRxnSubModel.genes);
[~, newRxnSubModel.rxnGeneMat, ~]=standardizeGrRules(newRxnSubModel);

%Since the newRxnSubModel is ready for incorportation, the iMK1208 model
%can be loaded for integration
load('iMK1208+suppInfo.mat','iMK1208');

%It should be noted that the incompatible nomenclatures (especially the
%metabolite identifiers) used in different GEMs and databases led to a
%serious problem in model comparison, curation and integration. In order to
%properly integrate the iMK1208 with this sub-model, both models should be
%unified into the same name space, since the metabolites in the sub-model
%use MetaCyc identifiers that are different from the ones used in iMK1208
%model. TableS2 contains the results for database mining and intensive
%manual curation while associating these metabolite identifiers. Load
%these manually selected reactions and their subSystems as the
%pre-processed array structure metaCycMetsIniMK.
load('iMK1208+suppInfo.mat','metaCycMetsIniMK');

%Replace the metabolites in the sub-model with the identifiers and names
%used in iMK1208 according to the mapping information in TableS2
%(metaCycMetsIniMK).
[a, b]=ismember(newRxnSubModel.mets,metaCycMetsIniMK);
I=find(a);
newRxnSubModel.mets=strcat(newRxnSubModel.mets,'_c');
newRxnSubModel.metNames(I)=iMK1208.metNames(b(I));
newRxnSubModel.mets(I)=iMK1208.mets(b(I));

%The necessary issues for merging iMK1208 with the new reactions determined
%from de novo generated draft GEMs have been resolved. Those can be
%directly merged into the Sco4 model, which refers to the 4th major update
%of the GEM for Streptomyces coelicolor.
Sco4=mergeModels({iMK1208 newRxnSubModel});

%Check out if the newly generated Sco4 could grow
sol=solveLP(Sco4);
disp(sol.f);
%The results indicate that this new model is functional
