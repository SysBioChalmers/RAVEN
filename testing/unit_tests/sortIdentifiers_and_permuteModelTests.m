%run this test case with the command
%results = runtests('sortIdentifiers_and_permuteModelTests.m')
function tests = sortIdentifiers_and_permuteModelTests
tests = functiontests(localfunctions);
end

function sortRandomizedModelTest(testCase)

% Load sorted toy model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat']);
taskStruct = parseTaskList(strcat(sourceDir, '/test_data/test_tasks.txt'));
rndModel = model;

% Randomly permutate model, do not use RAVEN functions
rndIdx = randperm(numel(rndModel.rxns));
fieldsToChange = {'rxns','lb','ub','rev','c','rxnNames','grRules','eccodes'};
for i=1:numel(fieldsToChange)
    rndModel.(fieldsToChange{i}) = rndModel.(fieldsToChange{i})(rndIdx);
end
rndModel.S          = rndModel.S(:,rndIdx);
rndModel.rxnGeneMat = rndModel.rxnGeneMat(rndIdx,:);

rndIdx = randperm(numel(rndModel.mets));
fieldsToChange = {'mets','metNames','metComps','metFormulas'};
for i=1:numel(fieldsToChange)
    rndModel.(fieldsToChange{i}) = rndModel.(fieldsToChange{i})(rndIdx);
end
rndModel.S     = rndModel.S(rndIdx,:);

rndIdx = randperm(numel(rndModel.genes));
fieldsToChange = {'genes','geneShortNames'};
for i=1:numel(fieldsToChange)
    rndModel.(fieldsToChange{i}) = rndModel.(fieldsToChange{i})(rndIdx);
end
rndModel.rxnGeneMat = rndModel.rxnGeneMat(:,rndIdx);

rndIdx = randperm(numel(rndModel.comps));
fieldsToChange = {'comps','compNames'};
for i=1:numel(fieldsToChange)
    rndModel.(fieldsToChange{i}) = rndModel.(fieldsToChange{i})(rndIdx);
end
[~,J]=sort(rndIdx);
[toreplace, bywhat] = ismember(rndModel.metComps,1:length(J));
rndModel.metComps(toreplace) = J(bywhat(toreplace));

%Sort random model
rndModel = sortIdentifiers(rndModel);

%Check that results is same as original model
verifyEqual(testCase,rndModel,model)
end
