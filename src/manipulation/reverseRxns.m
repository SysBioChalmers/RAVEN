function model=reverseRxns(model, rxns)
% reverseRxns
% Reverses the given reactions in place: negates their stoichiometry and
% objective coefficient, and swaps/negates their lower and upper bounds,
% so flux in the new positive direction is what was previously negative flux.
% model     the model to change (input/output)
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
