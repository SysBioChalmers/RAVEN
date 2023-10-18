%run this test case with the command
%results = runtests('modelSortingTests.m')
function tests = modelSortingTests
tests = functiontests(localfunctions);
end

function sortIdentifirs_and_permuteModelTest(testCase)

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
fieldsToChange = {'mets','metNames','metComps','metFormulas','metMiriams'};
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

function expandModel_and_contractModelTest(testCase)
%Load the expected model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

% Note that this does not work any model, as grRules might not be properly
% sorted (but should still be valid)
evalc('modelNew = expandModel(model);'); % Suppress warnings about complex grRules
modelNew = contractModel(modelNew);

verifyEqual(testCase,model,modelNew)
end

