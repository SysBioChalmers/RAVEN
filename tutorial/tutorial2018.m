%
% FILE NAME:    tutorial2018.m
% 
% DATE CREATED: 2018-01-08
%     MODIFIED: 2018-08-12
% 
% PROGRAMMER:   Hao Wang
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% 
% PURPOSE: This tutorial script demonstrates an automated process of
%          genome-scale model generation and curation for Streptomyces
%          coelicolor strin A3(2) by using the RAVEN Toolbox version 2.0.
%          This script includes two parts: de novo GEM reconstruciton and
%          model refinement based on manual curation.
% 
% 
% Input files:
%	  Sco_all_protein.faa   S. coelicolor proteome sequences in FASTA format
%           iMK1208.mat   The MATLAB file of iMK1208 model in RAVEN format
% SupportingTables.xlsx   Manual curation results in Excel format
%
% Output files:
% ScoCombinedModel.xlsx   The combined draft model in Excel format
%   newRxnSubModel.xlsx   The sub-model of new reacctions in Excel format
%

% Part One: De novo GEM Reconstruction

% 1. Draft GEM reconstruction using MetaCyc database

ScoMetaCycDraftModel=getMetaCycModelForOrganism('ScoMetaCyc','Sco_all_protein.faa',1);


% 2. Draft GEM reconstruction using KEGG database

% a. Draft GEM based on the genome information annotated by KEGG
ScoKEGGAnnotation=getKEGGModelForOrganism('sco','','','',0,0);
ScoKEGGAnnotation.id='ScoKEGGAnnotation';

% b. Draft GEM reconstructed from homology with HMMs traind from KEGG
% KO orthologs using the default cut-off values

% Define the type of  pre-trained hmms used in KEGG-based de novo reconstruction
dataDir = 'prok90_kegg82';

ScoKEGGHomology=getKEGGModelForOrganism('ScoKEGGHMMs',....
'Sco_all_protein.faa',dataDir,'',false,false);

% c. Direct merge the two KEGG-based  models
ScoKEGGDraftModel=mergeModels({ScoKEGGAnnotation ScoKEGGHomology});
% expandModel and contractModel for thorough integration
ScoKEGGDraftModel=expandModel(ScoKEGGDraftModel);
ScoKEGGDraftModel=contractModel(ScoKEGGDraftModel);
ScoKEGGDraftModel.id='mergedKEGG';
ScoKEGGDraftModel.description='Merged from reconstructions of KEGG annoation and homology search';


% 3. Combine the two draft GEMs obtained from KEGG and MetaCyc

ScoCombinedDraftModel=combineMetaCycKEGGModels(ScoMetaCycDraftModel, ScoKEGGDraftModel);

% Export the combined model into Excel format and subject for manual curation
exportToExcelFormat(ScoCombinedDraftModel,'ScoCombinedModel.xlsx');



% Part Two: Model Refinements Based on Manual Curation

% 4. A total of 398 reactions in the combined draft model were determined
% as new pathways that are absent from the iMK1208. They were extracted
% out as a sub model.

% a. Read in the list of manually selected reactions
[~, textData]=xlsread('SupportingTables.xlsx','TableS3');
selectedNewRxns.rxns=textData(2:end,1);
selectedNewRxns.subSystems=textData(2:end,3);
% b. Remove non-selected reactions from combined model and unused fields
rxnsToRemove=setdiff(ScoCombinedDraftModel.rxns,selectedNewRxns.rxns);
newRxnSubModel=removeReactions(ScoCombinedDraftModel,rxnsToRemove,true,true);
newRxnSubModel=rmfield(newRxnSubModel,{'geneFrom','grRulesKEGG','metFrom'});
% c. Update the subSystmes information
[a, b]=ismember(newRxnSubModel.rxns,selectedNewRxns.rxns);
newRxnSubModel.subSystems=selectedNewRxns.subSystems(b(find(a)));
% d. Define all new reactions in the Cytosol compartment 'c'
newRxnSubModel.comps{1}='c';
newRxnSubModel.compNames{1}='Cytoplasm';
newRxnSubModel.id='newRxnSubModel';
newRxnSubModel.description='The subModel of manually selected new reactions';
% e. Regenerate the genes field, because combined model has an empty rxnGeneMat
newRxnSubModel.genes={};
for k=1:numel(newRxnSubModel.rxns)
   newRxnSubModel.genes=[newRxnSubModel.genes;transpose(strsplit(newRxnSubModel.grRules{k},' or '))];
end
newRxnSubModel.genes=unique(newRxnSubModel.genes);
% f. Regenerate the rxnGeneMat field, and only 'or' relationship is considered
newRxnSubModel.rxnGeneMat=sparse(numel(newRxnSubModel.rxns),numel(newRxnSubModel.genes));
for i=1:numel(newRxnSubModel.rxns)
   %newRxnSubModel.genes=[newRxnSubModel.genes;transpose(strsplit(newRxnSubModel.grRules{k},' or '))];
   rxnGenes=strsplit(newRxnSubModel.grRules{i},' or ');
   [~, indexes]=ismember(rxnGenes,newRxnSubModel.genes);
   newRxnSubModel.rxnGeneMat(i,indexes)=1;
end
% g. Manually convert remaining KEGG mets to MetaCyc counterparts
% C00028 => Acceptor
metAIndex=find(strcmp(newRxnSubModel.mets,'C00028'));
metARxnIndex=find(newRxnSubModel.S(metAIndex,:));
metBIndex=find(strcmp(newRxnSubModel.mets,'Acceptor'));
newRxnSubModel.S(metBIndex,metARxnIndex)=newRxnSubModel.S(metAIndex,metARxnIndex);
% C00030 => Donor-H2
metAIndex=find(strcmp(newRxnSubModel.mets,'C00030'));
metARxnIndex=find(newRxnSubModel.S(metAIndex,:));
metBIndex=find(strcmp(newRxnSubModel.mets,'Donor-H2'));
newRxnSubModel.S(metBIndex,metARxnIndex)=newRxnSubModel.S(metAIndex,metARxnIndex);
% C03564 => DELTA1-PYRROLINE_2-CARBOXYLATE
metAIndex=find(strcmp(newRxnSubModel.mets,'C03564'));
metARxnIndex=find(newRxnSubModel.S(metAIndex,:));
metBIndex=find(strcmp(newRxnSubModel.mets,'DELTA1-PYRROLINE_2-CARBOXYLATE'));
newRxnSubModel.S(metBIndex,metARxnIndex)=newRxnSubModel.S(metAIndex,metARxnIndex);
newRxnSubModel=removeMets(newRxnSubModel,{'C00028','C00030','C03564'});

exportToExcelFormat(newRxnSubModel,'newRxnSubModel.xlsx');
%newRxnSubModel: #Rxns 398, #Mets 650(187), #Gene 404(106)


% 5. Generate the subModel of retrieved spontaneous reactions
% a. Load MetaCyc reactions and metabolites for model refinements
load('metaCycRxns.mat');
load('metaCycMets.mat');
% b. Read in the mapped MetaCyc reactions in iMK1208
[~, textData]=xlsread('SupportingTables.xlsx','TableS1');
metaCycRxnsIniMK=textData(3:end,8);
metaCycRxnsIniMK=metaCycRxnsIniMK(~cellfun(@isempty, metaCycRxnsIniMK));  %Remove empty elements
% c. Read in the mapped MetaCyc metabolites in iMK1208
[~, textData]=xlsread('SupportingTables.xlsx','TableS2');
metaCycMetsIniMK=textData(2:end,4);
keggMetsIniMK=textData(2:end,3);
metaCycMetsIniMK=metaCycMetsIniMK(~cellfun(@isempty, metaCycMetsIniMK));  %Remove empty elements
% d. Retrieve associated spontaneous reactions and obtain the submodel
% Get the list of all MetaCyc reactions
rxnList=unique([metaCycRxnsIniMK;newRxnSubModel.rxns]);
% Get the list of metabolites with MetaCyc ids
metList=unique([metaCycMetsIniMK;newRxnSubModel.mets]);
% Get the associated spontaneous reactions
spRxnList=retrieveSpontaneous(rxnList,metList);
% e. Generate spontaneous reaction submodel
rxnsToRemove=setdiff(metaCycRxns.rxns,spRxnList);
spRxnSubModel=removeReactions(metaCycRxns,rxnsToRemove,true,true);
% f. Regenerate the empty genes field
spRxnSubModel.genes{1}='s0001';
% g. Regenerate the rxnGeneMat with all zeros
spRxnSubModel.rxnGeneMat=ones(numel(spRxnSubModel.rxns),1);
spRxnSubModel.grRules=cell(numel(spRxnSubModel.rxns),1);
spRxnSubModel.grRules(:)={'s0001'};
spRxnSubModel.rxnConfidenceScores=cell(numel(spRxnSubModel.rxns),1);
spRxnSubModel.rxnConfidenceScores(:)={'2'};
spRxnSubModel.metNames=cell(numel(spRxnSubModel.mets),1);
spRxnSubModel.metCharges=zeros(numel(spRxnSubModel.mets),1);
[a, b]=ismember(spRxnSubModel.mets,metaCycMets.mets);
I=find(a);
spRxnSubModel.metNames(I)=metaCycMets.metNames(b(I));
spRxnSubModel.metCharges(I)=metaCycMets.metCharges(b(I));
I=cellfun(@isempty,spRxnSubModel.metNames);
spRxnSubModel.metNames(I)=spRxnSubModel.mets(I);
%spRxnSubModel: #Rxns 34(34), #Mets 77, #Gene 1


% 6. Prepare the metabolites for merge
% a. Load the iMK1208 model
model=load('iMK1208.mat');
% b. Match metabolites in newRxnSubModel to iMK1208
[a, b]=ismember(newRxnSubModel.mets,metaCycMetsIniMK);
I=find(a);
% Replace their metNames and ids with iMK ones
newRxnModel4Merge=newRxnSubModel;
for i=1:numel(newRxnModel4Merge.mets)
	newRxnModel4Merge.mets{i}=strcat(newRxnModel4Merge.mets{i},'_c');
end
newRxnModel4Merge.metNames(I)=model.metNames(b(I));
newRxnModel4Merge.mets(I)=model.mets(b(I));
% c. Match metabolites in spRxnSubModel to iMK
[hits, indexes]=ismember(spRxnSubModel.mets,metaCycMetsIniMK);
withits=find(hits);
% Replace their metNames and ids with iMK ones
spRxnModel4Merge=spRxnSubModel;
for i=1:numel(spRxnModel4Merge.mets)
	spRxnModel4Merge.mets{i}=strcat(spRxnModel4Merge.mets{i},'_c');
end
spRxnModel4Merge.metNames(withits)=model.metNames(indexes(withits));
spRxnModel4Merge.mets(withits)=model.mets(indexes(withits));
% d. Match metabolites in spRxnSubModel to newRxnSubModel
[hits, indexes]=ismember(spRxnModel4Merge.mets,newRxnModel4Merge.mets);
withits=find(hits);
spRxnModel4Merge.metNames(withits)=newRxnModel4Merge.metNames(indexes(withits));
spRxnModel4Merge.comps=newRxnModel4Merge.comps;
spRxnModel4Merge.compNames=newRxnModel4Merge.compNames;
spRxnModel4Merge.metComps=ones(numel(spRxnModel4Merge.mets),1);


% 7. Generate Sco4 model by adding up new metabolic and spontaneous reactions
Sco4=mergeModels({model newRxnModel4Merge spRxnModel4Merge});


% 8. Gap-filling for missing gene-associations in iMK1208
% a. New gene-associations from RAVEN2 reconstruction
[~, textData]=xlsread('SupportingTables.xlsx','TableS4');
gapfilling.rxns=textData(3:end,2);
gapfilling.grRules=textData(3:end,4);
% b. New gene-association for transport reactions
[~, textData]=xlsread('SupportingTables.xlsx','TableS6');
gapfilling.rxns=[gapfilling.rxns;textData(2:6,2)];
gapfilling.grRules=[gapfilling.grRules;textData(2:6,3)];
% c. Use new grRules from RAVEN2 reconstruction to fill gaps in iMK1208
[~, index]=ismember(gapfilling.rxns,Sco4.rxns);
Sco4.grRules(index)=gapfilling.grRules;
% d. Regenerate the rxnGeneMat
gapfilling.genes={};
for i=1:length(gapfilling.rxns)
	gapfilling.genes=[gapfilling.genes;transpose(strsplit(gapfilling.grRules{i},' or '))];
end
Sco4.genes=unique([Sco4.genes;gapfilling.genes]);
Sco4.geneMiriams=cell(numel(Sco4.genes),1);
Sco4.geneShortNames=Sco4.genes;
Sco4.rxnGeneMat=getRxnGeneMat(Sco4);


% 9. Add new transport reactions from MetaCyc reconstruction
% based on the curated transport reactions in Table S6
% a. First add new metabolites: gtca3[e], ENTEROBACTIN[c] and ENTEROBACTIN[e]
A=find(strcmp('gtca3_c',Sco4.mets));  % Find index for gtca3[c]
B=find(strcmp('ENTEROBACTIN',metaCycMets.mets)); % Find index for ENTEROBACTIN
metToAdd.mets={'gtca3_e';'ENTEROBACTIN_c';'ENTEROBACTIN_e'};
metToAdd.metNames{1,1}=Sco4.metNames{A};
metToAdd.metNames{2,1}=metaCycMets.metNames{B};
metToAdd.metNames{3,1}=metaCycMets.metNames{B};
metToAdd.compartments={'e';'c';'e'};
metToAdd.inchis{1,1}=Sco4.inchis{A};
metToAdd.inchis{2,1}=metaCycMets.inchis{B};
metToAdd.inchis{3,1}=metaCycMets.inchis{B};
metToAdd.metFormulas{1,1}=Sco4.metFormulas{A};
metToAdd.metFormulas{2,1}=metaCycMets.metFormulas{B};
metToAdd.metFormulas{3,1}=metaCycMets.metFormulas{B};
metToAdd.metMiriams{1,1}=Sco4.metMiriams{A};
metToAdd.metMiriams{2,1}=metaCycMets.metMiriams{B};
metToAdd.metMiriams{3,1}=metaCycMets.metMiriams{B};
metToAdd.metCharges(1,1)=Sco4.metCharges(A);
metToAdd.metCharges(2,1)=metaCycMets.metCharges(B);
metToAdd.metCharges(3,1)=metaCycMets.metCharges(B);
Sco4=addMets(Sco4,metToAdd,0);
% b. Prepare a structure of new transport reactions (Table S6)
rxnToAdd.rxns=textData(7:end,1);
rxnToAdd.grRules=textData(7:end,3);
rxnToAdd.subSystems=textData(7:end,4);
rxnToAdd.rxnNames=textData(7:end,5);
rxnToAdd.equations=textData(7:end,6);
rxnToAdd.eccodes=textData(7:end,8);
rxnToAdd.lb=[0; 0; 0; 0; 0];
rxnToAdd.ub=[1000; 1000; 1000; 1000; 1000];
rxnToAdd.confidenceScores={'2';'2';'2';'2';'2'};
% c. Add genes for these new transport reactions
geneToAdd.genes={};
for i=1:length(rxnToAdd.rxns)
	geneToAdd.genes=[geneToAdd.genes;transpose(strsplit(rxnToAdd.grRules{i},' or '))];
end
Sco4=addGenes(Sco4,geneToAdd);
% d. Add the new transport reactions
Sco4=addRxns(Sco4,rxnToAdd,1,'',0);


% 10. Add new exchange reactions
% Add sink reactions for the newly includeded secondary metabolites
Sco4=addExchangeRxns(Sco4,'out','CPD-14016_c');  % 2-methylisoborneol (RXN-12976)
Sco4=addExchangeRxns(Sco4,'out','CPD-9962_c');   % Albaflavenone (RXB-9352,RXN-9353)
Sco4=addExchangeRxns(Sco4,'out','CPD-11957_c');  % Desferrioxamine E (RXN-10998)
Sco4=addExchangeRxns(Sco4,'out','CPD-10680_c');  % 3,3-biflaviolin (RXN-9932)
Sco4=addExchangeRxns(Sco4,'out','CPD-10681_c');  % 3,8-biflaviolin (RXN-9933)
Sco4=addExchangeRxns(Sco4,'out','CPD-10158_c');  % Gemsmin (RXN-9476)
Sco4=addExchangeRxns(Sco4,'out','ACETONE_c');    % Gemsmin (RXN-9476)
Sco4=addExchangeRxns(Sco4,'out','CPD-14522_c');  % Hopanoid (RXN-13528)


% 11. Miscellaneous 
% a. Updated metabolite identifiers
newMets.metacyc=replace(Sco4.mets(1437:end),{'_c','_e'},'');
newMets.kegg=cell(numel(newMets.metacyc),1);
[a, b]=ismember(newMets.metacyc,metaCycMets.mets);
I=find(a);
newMets.kegg(I)=metaCycMets.keggid(b(I));
Sco4.metKEGGID=[model.metKEGGID;newMets.kegg];
Sco4.metMetaCycID=[metaCycMetsIniMK;newMets.metacyc];
% b. Reformat confidence scores
Sco4.rxnConfidenceScores(1860:end)={Sco4.rxnConfidenceScores{1859}};


% 12. Add model information
Sco4.id='Sco4';
Sco4.description='Streptomyces coelicolor genome-scale model';
% Annotation
annotation.defaultLB=-1000;
annotation.defaultUB=1000;
annotation.givenName='Hao';
annotation.familyName='Wang';
annotation.email='hao.wang@chalmers.se';
annotation.organization='Chalmers University of Technology';
annotation.taxonomy='100226';
Sco4.annotation=annotation;
