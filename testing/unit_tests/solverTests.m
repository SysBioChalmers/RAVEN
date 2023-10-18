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
    % Try all three types of flux minimization
    evalc('sol=solveLP(model,3);');    
    evalc('sol=solveLP(model,1);');
    evalc('sol=solveLP(model,0);');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    error('Solver not working')
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solOut');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solOut,'AbsTol',0.1) %Quite generous tolerance, as shadow price calculations fluctuate quite a bit between solvers
end

function testGurobi(testCase)
if exist('gurobi','file')~=3
    error('Gurobi not installed or cannot be found in MATLAB path, test skipped')
end
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('gurobi');

try
    % Try all three types of flux minimization
    evalc('sol=solveLP(model,3);');    
    evalc('sol=solveLP(model,1);');
    evalc('sol=solveLP(model,0);');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    error('Solver not working')
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solOut');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solOut,'AbsTol',0.1)
end

function testSoplex(testCase)
[a,~] = system('soplex --help');
if ~(a == 0)
    error('SoPlex not installed or cannot be found, test skipped')
end
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
setRavenSolver('soplex');

try
    % Try all three types of flux minimization
    evalc('sol=solveLP(model,3);');    
    evalc('sol=solveLP(model,1);');
    evalc('sol=solveLP(model,0);');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    error('Solver not working')
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solOutSoplex');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solOutSoplex,'AbsTol',0.1)
end

function testCobra(testCase)
if exist('initCobraToolbox.m','file')~=2
    error('COBRA Toolbox not installed or cannot be found in MATLAB path, test skipped')
end
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
try
    oldSolver=getpref('RAVEN','solver');
catch
end
global CBT_LP_SOLVER
global CBT_MILP_SOLVER
CBT_LP_SOLVER = 'glpk';
CBT_MILP_SOLVER = 'glpk';
setRavenSolver('cobra');

try
    % Try all three types of flux minimization
    evalc('sol=solveLP(model,3);');    
    evalc('sol=solveLP(model,1);');
    evalc('sol=solveLP(model,0);');
catch
    try
        setRavenSolver(oldSolver);
    catch
        rmpref('RAVEN','solver');
    end
    error('Solver not working')
end
try
    setRavenSolver(oldSolver);
catch
    rmpref('RAVEN','solver');
end

load([sourceDir,'/test_data/solverTestOutput.mat'], 'solOut');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solOut,'AbsTol',0.1)
end
