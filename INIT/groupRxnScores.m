function newRxnScores = groupRxnScores(model, origRxnScores, origRxnIds, groupIds, origRxnsToZero)
% groupRxnScores
% This function sums up the reaction scores for all reactions that were merged 
% into one by the linear merge.
% 
% model          The model with linearly merged rxns.
% origRxnScores  The rxnScores from the model before the linear merge.
% origRxnIds     The rxn ids of the model before the linear merge.
% groupIds       The groupIds vector output from linearMerge.
%                There is one integer for each rxn in origRxnIds. 0 means
%                the reaction was merged. A non-zero integer means that the
%                reaction was merged with all other rxns having the same integer.
% origRxnsToZero A logical vector saying which of the original rxns that should not
%                be part of the problem. The way this is solved is that all 
%                such reactions have a rxnScore of 0. If any original rxnScore 
%                value should be zero (which is very unlikely) it is changed to 0.01.
%                If the sum of the rxnScores for a merged rxn becomes zero 
%                while some of them are nonzero, the new value will also be 0.01, 
%                to distinguish the rxn from rxns with only rxns to zero.
%                There are two reasons why we don't want zeros in the reaction 
%                scores unless these reactions should be ignored:
%                  1) we want to be able to separate those
%                  2) it is difficult to handle a zero value in the MILP - 
%                     the on/off of such a reaction can be random, so better 
%                     to fix it in one direction.


newRxnScores = zeros(length(model.rxns),1);
[~,ia,ib] = intersect(model.rxns,origRxnIds);
grpIdsMerged = nan(length(model.rxns),1);
grpIdsMerged(ia) = groupIds(ib);
%check if any of the original scores are 0, in that case change them to 0.01 (unlikely)
origRxnScores(origRxnScores == 0) = 0.01;
%Then set the rxn scores for rxns to zero to 0
origRxnScores(origRxnsToZero) = 0;

%fill in original scores
newRxnScores(ia) = origRxnScores(ib);

for i = 1:length(model.rxns)
    %for reactions that are not merged with anything, just keep score as it is
    if grpIdsMerged(i) ~= 0
        %find all original rxns in the group
        sel = groupIds == grpIdsMerged(i);
        newRxnScores(i) = sum(origRxnScores(sel));
        if (newRxnScores(i) == 0 && any(origRxnScores(sel) ~= 0))
           %special unfortunate case, where the reactions happen to sum to 0 while some of them are nonzero
           %set to 0.01 in this case (proabably pretty unusual)
           newRxnScores(i) = 0.01;
        end
    end
end

end
