%run this test case with the command
%results = runtests('fillGapsSmallTests.m')
function tests = fillGapsSmallTests
tests = functiontests(localfunctions);
end

function testSmallSCIP(testCase)
%Test using small model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
modelDB=model; % Keep as database with reactions
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('scip');

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
    error('Solver not working')
end
sol=solveLP(model);
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end
%Should give non-zero flux
verifyTrue(testCase,-sol.f>0);
end

function testSmallGurobi(testCase)
if exist('gurobi','file')~=3
    error('Gurobi not installed or cannot be found in MATLAB path, test skipped')
end
%Test using small model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
modelDB=model; % Keep as database with reactions
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('gurobi');

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
    error('Solver not working')
end
sol=solveLP(model);
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end
%Expect at least 5% of the original growth
verifyTrue(testCase,-sol.f>0);
end