%run this test case with the command
%results = runtests('solverTests.m')
function tests = solverTests
tests = functiontests(localfunctions);
end

function testGlpk(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('glpk');

try
    sol=solveLP(model);
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    return
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solGlpk');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solGlpk,'AbsTol',1e-7)
end

function testGurobi(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('gurobi');

try
    sol=solveLP(model);
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    return
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solGurobi');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solGurobi,'AbsTol',1e-7)
end

function testCobra(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
global CBT_LP_SOLVER
CBT_LP_SOLVER = 'glpk';
setRavenSolver('cobra');

try
    sol=solveLP(model);
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    return
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solCobra');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solCobra,'AbsTol',1e-6)
end
