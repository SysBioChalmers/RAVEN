function allRxns=getAllRxnsFromGenes(model,rxns)
% getAllRxnsFromGenes
%   Given a list of reactions, this function finds the associated genes in
%   the template model and gives all reactions that are annotated by these
%   genes.
%
%   model       a model structure
%   rxns        either a cell array of IDs, a logical vector with the
%               same number of elements as reactions in the model, or a
%               vector of indexes
%
%   allRxns     either a cell array of IDs, a logical vector with the
%               same number of elements as reactions in the model, or a
%               vector of indexes, dependent on the format of rxns
%
%   Usage: allRxns=getAllRxnsFromGenes(model,rxns)

%If the supplied object is a character array, then convert it to a cell
%array
if ischar(rxns)
    rxns={rxns};
end

rxnIdx=getIndexes(model,rxns,'rxns');
[~, geneIdx]=find(model.rxnGeneMat(rxnIdx,:));
[allRxns, ~]=find(model.rxnGeneMat(:,geneIdx));
if iscell(rxns)
    allRxns=unique([model.rxns(allRxns);model.rxns(rxnIdx)]);
elseif isnumeric(rxns) || islogical(rxns)
    allRxns=unique([allRxns;rxnIdx]);
    if islogical(rxns)
        temp=false(numel(model.rxns),1);
        temp(allRxns)=true;
        allRxns=temp;
    end
end
end