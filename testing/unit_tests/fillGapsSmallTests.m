%run this test case with the command
%results = runtests('fillGapsSmallTests.m')
function tests = fillGapsSmallTests
tests = functiontests(localfunctions);
testGurobi = exist('gurobi','file')==3;
if testGurobi
    try
        gurobi_read('solverTests.m');
    catch ME
        if ~startsWith(ME.message,'Gurobi error 10012') % Expected error code, others may indicate problems with license
           testGurobi = false;
        end
    end
end
if ~testGurobi
    disp('Gurobi not installed or not functional, some fillGapsSmallTests skipped.')
    skipTests = contains({tests.Name},'gurobi','IgnoreCase',true);
    tests(skipTests) = [];
end
if exist('scip','file')~=3
    disp('SCIP MEX binary not installed or not functional, some fillGapsSmallTests skipped.')
    skipTests = contains({tests.Name},'scip','IgnoreCase',true);
    tests(skipTests) = [];
end
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
verifyTrue(testCase,sol.f>0);
end

function testSmallGurobi(testCase)
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
verifyTrue(testCase,sol.f>0);
end
