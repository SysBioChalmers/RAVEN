function [essentialRxns, essentialRxnsIndexes]=getEssentialRxns(model,ignoreRxns)
% getEssentialRxns
%   Calculate the essential reactions for a model to be solvable
%
%   model                   a model structure
%   ignoreRxns              cell array of reaction IDs which should not be
%                           checked (opt, default {})
%
%   essentialRxns           cell array with the IDs of the essential reactions
%   essentialRxnsIndexes    vector with the indexes of the essential reactions
%
%   Essential reactions are those which, when constrained to 0, result in an
%   infeasible problem.
%
%   Usage: [essentialRxns, essentialRxnsIndexes]=getEssentialRxns(model,ignoreRxns)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<2
    ignoreRxns={};
end

%Too make sure that it doesn't try to optimize for something
model.c=zeros(numel(model.rxns),1);

%First check that the problem is solvable
[sol, hsSolOut]=solveLP(model,1);

if sol.stat==-1 || isempty(sol.x)
    EM='No feasible solution to the full model';
    dispEM(EM);
end

%Check which reactions have flux. Only those can be essential. This is not
%the smallest list of reactions, but it's a fast way
rxnsToCheck=setdiff(model.rxns(abs(sol.x)>10^-8),ignoreRxns);
nToCheck=numel(rxnsToCheck);
minimize=true;
while 1
    if minimize==true
        sol=solveLP(setParam(model,'obj',rxnsToCheck,-1));
    else
        sol=solveLP(setParam(model,'obj',rxnsToCheck,1));
    end
    rxnsToCheck=intersect(rxnsToCheck,model.rxns(abs(sol.x)>10^-8));
    if numel(rxnsToCheck)>=nToCheck
        if minimize==true
            minimize=false;
        else
            break;
        end
    else
        nToCheck=numel(rxnsToCheck);
    end
end

essentialRxns={};
for i=1:numel(rxnsToCheck)
    sol=solveLP(setParam(model,'eq',rxnsToCheck(i),0),0,[],hsSolOut);
    if sol.stat==-1 || isempty(sol.x)
        essentialRxns=[essentialRxns;rxnsToCheck(i)];
    end
end

[~, essentialRxnsIndexes]=ismember(essentialRxns,model.rxns);
end
