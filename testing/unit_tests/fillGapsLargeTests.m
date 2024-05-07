%run this test case with the command
%results = runtests('fillGapsLargeTests.m')
function tests = fillGapsLargeTests
tests = functiontests(localfunctions);
testGurobi = exist('gurobi','file')==3;
if testGurobi
    try
        gurobi_read('test');
    catch ME
        testGurobi = false;
    end
end
if ~testGurobi
    disp('Gurobi not installed or cannot be found in MATLAB path, some fillGapsLargeTests skipped.')
    skipTests = contains({tests.Name},'gurobi','IgnoreCase',true);
    tests(skipTests) = [];
end
if exist('scip','file')~=3
    disp('SCIP MEX binary not installed or not functional, test skipped.')
    skipTests = contains({tests.Name},'scip','IgnoreCase',true);
    tests(skipTests) = [];
end
end

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

%Remove first 10 reactions
model=removeReactions(modelDB,(1:10));
modelDB.id='DB';
try
    evalc('[newConnected,cannotConnect,addedRxns,model,exitFlag]=fillGaps(model,modelDB,false,false)');
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

function testLargeSCIP(testCase)
sourceDir = fileparts(fileparts(fileparts(which(mfilename))));
evalc('model=importModel(fullfile(sourceDir,''tutorial'',''iAL1006 v1.00.xml''))');
model.c(1484)=1;
modelDB=model; % Keep as database with reactions
% Force growth in gapped model
sol=solveLP(model);
model.lb(1484)=abs(sol.f*0.1);
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('scip');

%Remove first 10 reactions
model=removeReactions(model,(1:10));
modelDB.id='DB';
try
    evalc('[newConnected,cannotConnect,addedRxns,model,exitFlag]=fillGaps(model,modelDB,false,true)');
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
