%run this test case with the command
%results = runtests('sortIdentifiers_and_permuteModelTests.m')
function tests = sortIdentifiers_and_permuteModelTests
tests = functiontests(localfunctions);
end

function sortRandomizedModelTest(testCase)

%Load the expected (i.e. sorted) model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
expModel = model;

%Create the actual model that will be permuted and sorted
actModel = expModel;

%Randomly permutate model, do not use RAVEN functions
rndIdx = randperm(numel(actModel.rxns));
fieldsToChange = {'rxns','lb','ub','rev','c','rxnNames','grRules','eccodes'};
for i=1:numel(fieldsToChange)
    actModel.(fieldsToChange{i}) = actModel.(fieldsToChange{i})(rndIdx);
end
actModel.S          = actModel.S(:,rndIdx);
actModel.rxnGeneMat = actModel.rxnGeneMat(rndIdx,:);

rndIdx = randperm(numel(actModel.mets));
fieldsToChange = {'mets','metNames','metComps','metFormulas'};
for i=1:numel(fieldsToChange)
    actModel.(fieldsToChange{i}) = actModel.(fieldsToChange{i})(rndIdx);
end
actModel.S     = actModel.S(rndIdx,:);

rndIdx = randperm(numel(actModel.genes));
fieldsToChange = {'genes','geneShortNames'};
for i=1:numel(fieldsToChange)
    actModel.(fieldsToChange{i}) = actModel.(fieldsToChange{i})(rndIdx);
end
actModel.rxnGeneMat = actModel.rxnGeneMat(:,rndIdx);

rndIdx = randperm(numel(actModel.comps));
fieldsToChange = {'comps','compNames'};
for i=1:numel(fieldsToChange)
    actModel.(fieldsToChange{i}) = actModel.(fieldsToChange{i})(rndIdx);
end
[~,J]=sort(rndIdx);
[toreplace, bywhat] = ismember(actModel.metComps,1:length(J));
actModel.metComps(toreplace) = J(bywhat(toreplace));

%Sort randomly permutated model
actModel = sortIdentifiers(actModel);

%Check that the actual model is the same as the expected model
verifyEqual(testCase,actModel,expModel)
end
