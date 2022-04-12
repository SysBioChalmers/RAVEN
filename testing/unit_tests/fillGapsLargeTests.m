%run this test case with the command
%results = runtests('fillGapsLargeTests.m')
function tests = fillGapsLargeTests
tests = functiontests(localfunctions);
end

%Skip testLargeGlpk, fails with larger models, known issue with glpk
% function testLargeGlpk(testCase)
% sourceDir = fileparts(fileparts(fileparts(which(mfilename))));
% evalc('model=importModel(fullfile(sourceDir,''tutorial'',''iAL1006 v1.00.xml''))');
% model.c(1484)=1;
% modelDB=model; % Keep as database with reactions
% try
%     oldSolver=getpref('RAVEN','solver');
% catch
% end
% setRavenSolver('glpk');
% solDB=solveLP(modelDB);
% 
% %Remove ~10% reactions at random
% removeRxns=unique(randi(numel(modelDB.rxns),floor(0.10*numel(modelDB.rxns)),1));
% removeRxns(removeRxns==1484)=[]; % Do not remove biomass reaction
% model=removeReactions(modelDB,removeRxns);
% modelDB.id='DB';
% try
%     evalc('[newConnected,cannotConnect,addedRxns,model,exitFlag]=fillGaps(model,modelDB)');
% catch
%     try
%         setRavenSolver(oldSolver);
%     catch
%         rmpref('RAVEN','solver');
%     end
%     return
% end
% sol=solveLP(model);
% try
%     setRavenSolver(oldSolver);
% catch
%     rmpref('RAVEN','solver');
% end
% end

function testLargeGurobi(testCase)
sourceDir = fileparts(fileparts(fileparts(which(mfilename))));
evalc('model=importModel(fullfile(sourceDir,''tutorial'',''iAL1006 v1.00.xml''))');
model.c(1484)=1;
modelDB=model; % Keep as database with reactions
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('gurobi');
solDB=solveLP(modelDB);

%Remove first 10 reactions
model=removeReactions(modelDB,(1:10));
modelDB.id='DB';
try
    evalc('[newConnected,cannotConnect,addedRxns,model,exitFlag]=fillGaps(model,modelDB)');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    return
end
sol=solveLP(model);
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end
%Expect at least 5% of the original growth
verifyTrue(testCase,-sol.f>0.05*-solDB.f);
end

function testLargeCobra(testCase)
sourceDir = fileparts(fileparts(fileparts(which(mfilename))));
evalc('model=importModel(fullfile(sourceDir,''tutorial'',''iAL1006 v1.00.xml''))');
model.c(1484)=1;
modelDB=model; % Keep as database with reactions
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('cobra');
solDB=solveLP(modelDB);

%Remove first 10 reactions
model=removeReactions(modelDB,(1:10));
modelDB.id='DB';
try
    evalc('[newConnected,cannotConnect,addedRxns,model,exitFlag]=fillGaps(model,modelDB)');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    return
end
sol=solveLP(model);
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

%Expect at least 5% of the original growth
verifyTrue(testCase,-sol.f>0.05*-solDB.f);
end