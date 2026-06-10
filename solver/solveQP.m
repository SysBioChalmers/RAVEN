function solution=solveQP(model,rxns,values,maxIter, restartIter)
% solveQP  Solve a quadratic fitting problem.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxns : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of indexes
%     to fit to.
% values : double
%     the values to fit the fluxes to.
% maxIter : double, optional
%     maximal number of iterations (default 1000).
% restartIter : double, optional
%     run the fitting up to this many times in case it does not converge
%     (default 1).
%
% Returns
% -------
% solution : struct
%     solution with fields:
%
%     - f : objective value.
%     - x : primal.
%     - stat : exit flag.
%
% Examples
% --------
%     solution = solveQP(model, rxns, values, maxIter, restartIter);
rxns=char(rxns);

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
