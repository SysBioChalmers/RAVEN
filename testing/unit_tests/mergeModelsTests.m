%run this test case with the command
%results = runtests('mergeModelsTests.m')
function tests = mergeModelsTests
tests = functiontests(localfunctions);
end

function mergeModelsBasicTest(testCase)
%Merging plain models (without rxnFrom/metFrom/geneFrom) should create these
%fields and trace each reaction/metabolite/gene back to the model it came from
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelA    = model;
modelA.id = 'modelA';

%Make all reactions and metabolites in modelB unique so that they are all
%added as new entities when merging (metabolites are matched by name)
modelB          = model;
modelB.id       = 'modelB';
modelB.rxns     = strcat(modelB.rxns,'_B');
modelB.mets     = strcat(modelB.mets,'_B');
modelB.metNames = strcat(modelB.metNames,'_B');

evalc('modelMerged=mergeModels({modelA;modelB})');

%The From fields should exist and have the same size as their entities. Before
%the fix, metFrom was not extended and ended up shorter than mets
verifyTrue(testCase,isfield(modelMerged,'rxnFrom'));
verifyTrue(testCase,isfield(modelMerged,'metFrom'));
verifyTrue(testCase,isfield(modelMerged,'geneFrom'));
verifyEqual(testCase,numel(modelMerged.rxnFrom),numel(modelMerged.rxns));
verifyEqual(testCase,numel(modelMerged.metFrom),numel(modelMerged.mets));
verifyEqual(testCase,numel(modelMerged.geneFrom),numel(modelMerged.genes));

%Each reaction and metabolite should be traced back to its source model
verifyEqual(testCase,modelMerged.rxnFrom,[repmat({'modelA'},numel(modelA.rxns),1);repmat({'modelB'},numel(modelB.rxns),1)]);
verifyEqual(testCase,modelMerged.metFrom,[repmat({'modelA'},numel(modelA.mets),1);repmat({'modelB'},numel(modelB.mets),1)]);
%Genes are shared between the models, so all are traced to the first model
verifyTrue(testCase,all(ismember(modelMerged.geneFrom,{'modelA','modelB'})));
end

function mergeModelsKeepExistingFromTest(testCase)
%If the input models already have rxnFrom/metFrom/geneFrom fields, mergeModels
%should keep those values instead of overwriting them with the model ids
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelA          = model;
modelA.id       = 'modelA';
modelA.rxnFrom  = repmat({'origA'},numel(modelA.rxns),1);
modelA.metFrom  = repmat({'origA'},numel(modelA.mets),1);

modelB          = model;
modelB.id       = 'modelB';
modelB.rxns     = strcat(modelB.rxns,'_B');
modelB.mets     = strcat(modelB.mets,'_B');
modelB.metNames = strcat(modelB.metNames,'_B');
modelB.rxnFrom  = repmat({'origB'},numel(modelB.rxns),1);
modelB.metFrom  = repmat({'origB'},numel(modelB.mets),1);

evalc('modelMerged=mergeModels({modelA;modelB})');

%The pre-existing labels should be preserved, not replaced by the model ids
verifyEqual(testCase,modelMerged.rxnFrom,[repmat({'origA'},numel(modelA.rxns),1);repmat({'origB'},numel(modelB.rxns),1)]);
verifyEqual(testCase,modelMerged.metFrom,[repmat({'origA'},numel(modelA.mets),1);repmat({'origB'},numel(modelB.mets),1)]);
end

function mergeModelsViaCopyToCompsTest(testCase)
%copyToComps runs mergeModels with the copyToComps flag set. When the input
%model has no rxnFrom/metFrom/geneFrom fields, the output should not gain any
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

evalc('modelNew=copyToComps(model,{''p''},''ACKr'')');

verifyFalse(testCase,isfield(modelNew,'rxnFrom'));
verifyFalse(testCase,isfield(modelNew,'metFrom'));
verifyFalse(testCase,isfield(modelNew,'geneFrom'));
%The reaction should still have been copied to the new compartment
verifyTrue(testCase,ismember('ACKr_p',modelNew.rxns));
end

function mergeModelsViaCopyToCompsKeepFromTest(testCase)
%When the input model already has the From fields, copyToComps should keep
%them and assign empty labels ('') to the newly created reactions/metabolites
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

model.rxnFrom  = repmat({'orig'},numel(model.rxns),1);
model.metFrom  = repmat({'orig'},numel(model.mets),1);
model.geneFrom = repmat({'orig'},numel(model.genes),1);

nRxns = numel(model.rxns);
nMets = numel(model.mets);

evalc('modelNew=copyToComps(model,{''p''},''ACKr'')');

%The fields should still be present and consistent in size
verifyEqual(testCase,numel(modelNew.rxnFrom),numel(modelNew.rxns));
verifyEqual(testCase,numel(modelNew.metFrom),numel(modelNew.mets));
%The original labels are kept, the new reactions/metabolites are labelled ''
verifyEqual(testCase,modelNew.rxnFrom(1:nRxns),model.rxnFrom);
verifyEqual(testCase,modelNew.metFrom(1:nMets),model.metFrom);
verifyTrue(testCase,all(strcmp(modelNew.rxnFrom(nRxns+1:end),'')));
verifyTrue(testCase,all(strcmp(modelNew.metFrom(nMets+1:end),'')));
%No new genes are introduced, so geneFrom is left untouched
verifyEqual(testCase,modelNew.geneFrom,model.geneFrom);
end

