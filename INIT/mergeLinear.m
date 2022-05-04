function [reducedModel,origRxnIds,groupIds,reversedRxns]=mergeLinear(model, noMergeRxns)
% mergeLinear
%   Simplifies a model by merging rxns with linear dependencies. There are two
%   differences from the function in simplifyModel:
%   1) Here we return which reactions were merged
%   2) Here we allow reversible reactions to be merged as well, which will 
%      probably reduce the size of the model more.
%
%   model                 a model structure
%   noMergeRxns           Cell array with reaction IDs that are not allowed to be
%                         merged
%
%   reducedModel          an updated model structure
%   origRxnIds            Cell array: original rxn ids, used together with the groupIds variable
%   groupIds              Vector of ids: merged rxns have the same id, same length as origRxnIds
%                         0 means a reaction was not merged.
%   reversedRxns          Vector of booleans, saying for each rxn if it is reversed in the combined rxns.
%
%   Usage: [reducedModel,origRxnIds,groupIds]=mergeLinear(model,noMergeRxns)


reducedModel=model;


reducedModel.genes={};
reducedModel.rxnGeneMat=sparse(numel(reducedModel.rxns),0);
reducedModel.grRules(:)={''};

if isfield(reducedModel,'geneShortNames')
    reducedModel.geneShortNames={};
end
if isfield(reducedModel,'geneMiriams')
    reducedModel.geneMiriams={};
end
if isfield(reducedModel,'geneComps')
    reducedModel.geneComps=[];
end

nextGroupId = 1;
origRxnIds = reducedModel.rxns;
groupIds = zeros(numel(reducedModel.rxns),1);
reversedRxns = false(numel(reducedModel.rxns),1);

%Loop through and iteratively group linear reactions
while 1
    %Get the banned reaction indexes. Note that the indexes will change
    %in each iteration, but the names will not as they won't be merged
    %with any other reaction
    bannedIndexes=getIndexes(reducedModel,noMergeRxns,'rxns');

    %Select all metabolites that are only present as reactants/products
    %in one reaction
    twoNonZero = find(sum(reducedModel.S ~= 0, 2) == 2);

    mergedSome=false;

    %Loop through each of them and see if the reactions should be
    %merged
    for i=1:numel(twoNonZero)
        involvedRxns=find(reducedModel.S(twoNonZero(i),:));
        %Check that we can have one positive and one negative
        pos = sum(reducedModel.S(twoNonZero(i),involvedRxns).' > 0 | reducedModel.rev(involvedRxns));
        neg = sum(reducedModel.S(twoNonZero(i),involvedRxns).' < 0 | reducedModel.rev(involvedRxns));
        
                
        %Check so that one or both of the reactions haven't been merged
        %already
        if numel(involvedRxns)==2 && isempty(intersect(bannedIndexes,involvedRxns)) && pos >= 1 && neg >= 1
            %first, take care of a special case: If the first reaction is producing the metabolite and if it is reversible, 
            %and the second is also producing it and is not reversible, change the order - the code below will not work otherwise
            if reducedModel.rev(involvedRxns(1)) && (~reducedModel.rev(involvedRxns(2))) && ...
               (reducedModel.S(twoNonZero(i),involvedRxns(1)) > 0) && (reducedModel.S(twoNonZero(i),involvedRxns(2)) > 0)
                involvedRxns = flip(involvedRxns);
            end
            
            %first make sure the first reaction is producing the metabolite
            if reducedModel.S(twoNonZero(i),involvedRxns(1)) < 0
                %it is not producing the metabolite - fix that
                %first choice: use the second reaction as producer if it is producing
                if reducedModel.S(twoNonZero(i),involvedRxns(2)) > 0
                    involvedRxns = flip(involvedRxns);%make the second the first
                else
                    %now we know that the second reaction is not producing, so we can safely try to make the first a producer
                    if reducedModel.rev(involvedRxns(1)) == 1
                        [reducedModel,reversedRxns] = flipRxn(reducedModel, involvedRxns(1), groupIds, reversedRxns);
                    else %ok, finally try to flip the second reaction
                        if reducedModel.rev(involvedRxns(2)) == 1
                            [reducedModel,reversedRxns] = flipRxn(reducedModel, involvedRxns(2), groupIds, reversedRxns);
                            involvedRxns = flip(involvedRxns);%make the second the first
                        else
                            error('We should never end up here!');
                        end
                    end
                end
            end
            %Now, make sure the second rxn is a consumer
            if reducedModel.S(twoNonZero(i),involvedRxns(2)) > 0
                if reducedModel.rev(involvedRxns(2)) == 1
                    [reducedModel,reversedRxns] = flipRxn(reducedModel, involvedRxns(2), groupIds, reversedRxns);
                else
                    error('We should never end up here!');
                end
            end

            %Calculate how many times the second reaction has to be
            %multiplied before being merged with the first
            stoichRatio=abs(reducedModel.S(twoNonZero(i),involvedRxns(1))/reducedModel.S(twoNonZero(i),involvedRxns(2)));

            %Add the second to the first
            reducedModel.S(:,involvedRxns(1))=reducedModel.S(:,involvedRxns(1))+reducedModel.S(:,involvedRxns(2))*stoichRatio;

            %Clear the second reaction
            reducedModel.S(:,involvedRxns(2))=0;

            %This is to prevent numerical issues. It should be 0
            %already
            reducedModel.S(twoNonZero(i),involvedRxns(1))=0;

            %At this point the second reaction is certain to be deleted
            %in a later step and can therefore be ignored

            %Recalculate the bounds for the new reaction. This can be
            %problematic since the scale of the bounds may change
            %dramatically. Let the most constraining reaction determine
            %the new bound
            lb1=reducedModel.lb(involvedRxns(1));
            lb2=reducedModel.lb(involvedRxns(2));
            ub1=reducedModel.ub(involvedRxns(1));
            ub2=reducedModel.ub(involvedRxns(2));

            if lb2~=-inf
                reducedModel.lb(involvedRxns(1))=max(lb1,lb2/stoichRatio);
            end
            if ub2~=inf
                reducedModel.ub(involvedRxns(1))=min(ub1,ub2/stoichRatio);
            end
            
            %take care of the .rev flag - it could be that the combined rxn changes from rev to irrev
            reducedModel.rev(involvedRxns(1)) = reducedModel.rev(involvedRxns(1))*reducedModel.rev(involvedRxns(2));%this is a way to do an "and" operation with 0 and 1 numbers

            %Then recalculate the objective coefficient. The resulting
            %coefficient is the weighted sum of the previous
            reducedModel.c(involvedRxns(1))=reducedModel.c(involvedRxns(1))+reducedModel.c(involvedRxns(2))*stoichRatio;

            %store which reactions that have been merged
            rxnInd1 = find(strcmp(origRxnIds, reducedModel.rxns(involvedRxns(1))));
            rxnInd2 = find(strcmp(origRxnIds, reducedModel.rxns(involvedRxns(2))));
            grpId = max(groupIds(rxnInd1),groupIds(rxnInd2));
            if grpId == 0
               grpId = nextGroupId;
               nextGroupId = nextGroupId + 1;
            end

            if groupIds(rxnInd1) ~= grpId
               if groupIds(rxnInd1) == 0
                   %not merged before, just set the group id
                   groupIds(rxnInd1) = grpId;
               else
                   %merged before - all rxns with the same group id should be changed
                   groupIds(groupIds == groupIds(rxnInd1)) = grpId;
               end
            end
            if groupIds(rxnInd2) ~= grpId
               if groupIds(rxnInd2) == 0
                   %not merged before, just set the group id
                   groupIds(rxnInd2) = grpId;
               else
                   %merged before - all rxns with the same group id should be changed
                   groupIds(groupIds == groupIds(rxnInd2)) = grpId;
               end
            end

            %Iterate again
            mergedSome=true;
        end
    end

    %All possible reactions merged
    if mergedSome==false
        break;
    end

    %Now delete all reactions that involve no metabolites
    I=find(sum(reducedModel.S~=0,1)==0);

    %Remove reactions
    reducedModel=removeReactions(reducedModel,I);

    %Remove metabolites
    notInUse=sum(reducedModel.S~=0,2)==0;
    reducedModel=removeMets(reducedModel,notInUse);
end
    
function [model1,reversedRxns1] = flipRxn(model1, rxnInd, groupIds1, reversedRxns1)
    model1.S(:,rxnInd) = model1.S(:,rxnInd)*-1;
    %swap the bounds
    ub = model1.ub(rxnInd);
    model1.ub(rxnInd) = -model1.lb(rxnInd);
    model1.lb(rxnInd) = -ub;
    %flip the objective
    model1.c(rxnInd) = -model1.c(rxnInd);
    
    %now take care of the reversedRxns - if this is a group, reverse all of the
    %reactions in the group in the reversedRxns index - they will all be reversed at the
    %same time since they are the same rxn.
    rxnIndices = rxnInd;
    if groupIds1(rxnInd) > 0
        rxnIndices = find(groupIds1 == groupIds1(rxnInd));
    end
    reversedRxns1(rxnIndices) = ~reversedRxns1(rxnIndices);
end
end
