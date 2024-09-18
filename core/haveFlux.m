function I=haveFlux(model,cutOff,rxns)
% haveFlux
%   Checks which reactions can carry a (positive or negative) flux. Is used
%   as a faster version of getAllowedBounds if it is only interesting
%   whether the reactions can carry a flux or not
%
% Input:
%   model       a model structure
%   cutOff      the flux value that a reaction has to carry to be
%               identified as positive (optional, default 10^-8)
%   rxns        either a cell array of IDs, a logical vector with the
%               same number of elements as metabolites in the model, or a
%               vector of indexes (optional, default model.rxns)
%
% Output:
%   I           logical array with true if the corresponding reaction can
%               carry a flux
%
% If a model has +/- Inf bounds then those are replaced with an arbitary
% large value of +/- 10000 prior to solving
%
% Usage: I = haveFlux(model, cutOff, rxns)

if nargin<2
    cutOff=10^-6;
end
if isempty(cutOff)
    cutOff=10^-6;
end
if nargin<3
    rxns=model.rxns;
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
end

%This is since we're maximizing for the sum of fluxes, which isn't possible
%when there are infinite bounds
model.lb(model.lb==-inf)=-10000;
model.ub(model.ub==inf)=10000;

%Get the reaction IDs. A bit of an awkward way, but fine.
indexes=getIndexes(model,rxns,'rxns');
rxns=model.rxns(indexes);

%First remove all dead-end reactions
smallModel=simplifyModel(model,false,false,true,true);

%Get the indexes of the reactions in the connected model
indexes=getIndexes(smallModel,intersect(smallModel.rxns,rxns),'rxns');
J=false(numel(indexes),1);
mixIndexes=indexes(randperm(numel(indexes)));

%Maximize for all fluxes first in order to get fewer rxns to test
smallModel.c=ones(numel(smallModel.c),1);
sol=solveLP(smallModel);
if ~isempty(sol.x)
    J(abs(sol.x(mixIndexes))>cutOff)=true;
end

%Loop through and maximize then minimize each rxn if it does not already
%have a flux
Z=zeros(numel(smallModel.c),1);
hsSolOut=[];
for i=[1 -1]
    for j=1:numel(J)
        if J(j)==false
            %Only minimize for reversible reactions
            if i==1 || smallModel.rev(mixIndexes(j))~=0
                smallModel.c=Z;
                smallModel.c(mixIndexes(j))=i;
                [sol, hsSolOut]=solveLP(smallModel,0,[],hsSolOut);
                if any(sol.x)
                    J(abs(sol.x(mixIndexes))>cutOff)=true;
                end
            end
        end
    end
end

I=ismember(rxns,smallModel.rxns(mixIndexes(J)));
end
