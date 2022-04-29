function model=reverseRxns(model, rxns)
% model=reverseRxns(model, rxns)
% reverses reactions
% model     the model to change
% rxns      the rxns to reverse

rxnInd = find(ismember(model.rxns,rxns));

%reverse the reactions
model.S(:,rxnInd) = model.S(:,rxnInd)*-1;
%swap the bounds
ub = model.ub(rxnInd);
model.ub(rxnInd) = -model.lb(rxnInd);
model.lb(rxnInd) = -ub;
%flip the objective
model.c(rxnInd) = -model.c(rxnInd);

end
