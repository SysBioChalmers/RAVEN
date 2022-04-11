%run this test case with the command
%results = runtests('solverTests.m')
function tests = solverTests
tests = functiontests(localfunctions);
end

function testGlpk(testCase)
%This function tests the glpk solver
%Load the expected (i.e. sorted) model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
oldSolver=getpref('RAVEN','solver');
setRavenSolver('glpk');
sol=solveLP(model);
setRavenSolver(oldSolver);
load([sourceDir,'/test_data/solverTestOutput.mat'], 'solGlpk');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solGlpk,'AbsTol',1e-7)
end

function testGurobi(testCase)
%This function tests the glpk solver
%Load the expected (i.e. sorted) model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
oldSolver=getpref('RAVEN','solver');
setRavenSolver('gurobi');
sol=solveLP(model);
setRavenSolver(oldSolver);
load([sourceDir,'/test_data/solverTestOutput.mat'], 'solGurobi');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solGurobi,'AbsTol',1e-7)
end

function testCobra(testCase)
%This function tests the glpk solver
%Load the expected (i.e. sorted) model
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
oldSolver=getpref('RAVEN','solver');
global CBT_LP_SOLVER
CBT_LP_SOLVER = 'glpk';
setRavenSolver('cobra');
sol=solveLP(model);
setRavenSolver(oldSolver);
load([sourceDir,'/test_data/solverTestOutput.mat'], 'solCobra');
%Check that the actual model is the same as the expected model
verifyEqual(testCase,sol,solCobra,'AbsTol',1e-6)
end
