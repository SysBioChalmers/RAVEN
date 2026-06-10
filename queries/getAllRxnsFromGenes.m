function allRxns=getAllRxnsFromGenes(model,rxns)
% getAllRxnsFromGenes  Find all reactions annotated by the genes of a set.
%
% Given a list of reactions, this function finds the associated genes in
% the model and returns all reactions that are annotated by these genes.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxns : cell or logical or double
%     either a cell array of IDs, a logical vector with the same number of
%     elements as reactions in the model, or a vector of indexes.
%
% Returns
% -------
% allRxns : cell or logical or double
%     either a cell array of IDs, a logical vector with the same number of
%     elements as reactions in the model, or a vector of indexes,
%     dependent on the format of rxns.
%
% Examples
% --------
%     allRxns = getAllRxnsFromGenes(model, rxns);

if ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
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