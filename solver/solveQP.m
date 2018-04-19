function solution=solveQP(model,rxns,values,maxIter, restartIter)
% solveQP
%   Solves a quadratic fitting problem.
%
%   model         a model structure
%   rxns          either a cell array of reaction IDs, a logical vector 
%                 with the same number of elements as reactions in the model,
%                 of a vector of indexes to fit to
%   values        the values to fit the fluxes to
%   maxIter       maximal number of iterations (opt, default 1000)
%   restartIter   run the fitting up to this many times in case it does
%                 not converge (opt, default 1)
%
%   solution
%         f       Objective value
%         x       Primal
%         stat    Exit flag
%
%   Usage: solution=solveQP(model,rxns,values,maxIter, restartIter)
%
%   Rasmus Agren, 2013-02-13
%

if nargin<4
    maxIter=1000;
end

if nargin<5
    restartIter=1;
end

%Check that it's feasible
solution=solveLP(model);
if isempty(solution.x)
    return;
end

%Get the indexes of the fluxes to fit
allIndexes=getIndexes(model,rxns,'rxns');

f=zeros(numel(model.rxns),1);
H=zeros(numel(model.rxns));

%Get the fitIndexes for the experiment
for j=1:numel(allIndexes) %Not that neat
    H(allIndexes(j),allIndexes(j))=2;
end

f(allIndexes)=values.*-2;

%For a badly formulated problem it can occur that the solver stalls. This
%can sometimes be fixed by trying to solve the problem again
options=optimset('MaxIter',maxIter);
for j=1:restartIter
    [x,fval,flag] = ...
        quadprog(H,f,[],[],model.S,model.b,model.lb,model.ub,[],options);
    if flag>0
        break;
    end
end

solution.f=fval;
solution.x=x;
solution.stat=flag;
end
