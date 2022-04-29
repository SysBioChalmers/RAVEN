function model=rescaleModelForINIT(model, maxStoichVal)
% model=rescaleModelForINIT(model, maxStoichVal)
%
% The idea with this function is to rescale the MILP problem in ftINIT to avoid large differences
% in flux magnitudes between reactions.
% For now it just scales down all reactions with high stoichiometric coefficients
% There is room for improvement here - the best would be to convert mets such as albumin
% to instead represent 1/100 albumin - that would create much less extreme coefficients.
%
% model         the model to be modified
% maxStoichVal  all reactions with stoichiometric coefficent higher than this 
%               will be scaled down. (opt, default 250)

if (nargin < 2)
    maxStoichVal = 250;
end


%find all reactions containing high stoichiometry coefficients
maxCoeffs = max(model.S, [], 1);
rxnsToCheck = maxCoeffs > maxStoichVal;

%constructEquations(model, model.rxns(rxnsToCheck))


%pointlessMets = {'H+';'H2O';'Pi'}

rxnsToCheckInd = find(rxnsToCheck);

%We just scale all reactions so they don't have a coefficient > 250 for now

%do this with a loop, doesn't really matter
for rxnInd = rxnsToCheckInd
    %find the lowest index which is not a pointless metabolite, such as H+, H2O etc.
    rxnS = model.S(:,rxnInd);
    %metInd = find(rxnS ~= 0);
    
    largestCoeff = max(abs(rxnS));
    scaleFactor = maxStoichVal/largestCoeff;
    model.S(:,rxnInd) = model.S(:,rxnInd) .* scaleFactor;
    
    
    %model.metNames(metInd)
    %filter out uninteresting mets
    %iToRem = [];
    %for i = 1:length(metInd)
    %   if ismember(model.metNames(metInd(i)), pointlessMets) 
    %       iToRem = [iToRem;i];
    %   end
    %end
    %metIndFilt = metInd;
    %metIndFilt(iToRem) = [];
    %[val,ind] = min(abs(rxnS(metIndFilt)))
end



end
%{
function [model,success] = redMet(model, rxnInd, metInd, scaleFactor)
    %The problem here is that some metabolites, such as large proteins,
    %are large - the reactions that build/degrade them will thus have very
    %high coefficients on metabolites like ATP and amino acids.
    %What we want to do is to change these reactions and replace the large metabolite
    %with a fraction of that metabolite, for example 0.001 albumin or something.
    %We only do this if the reactions involved are simple, i.e. transport reactions in and out
    %regarding the albumin. 
    
    allRxnsWithMet = find(model.S(metInd,:) ~=0);

    for i = allRxnsWithMet
        
    end
    
end

function b = checkMetRec(model, checkedRxns, metInd)
    b = true;
    allRxnsWithMet = find(model.S(metInd,:) ~=0);
    checkedRxnsNextLevel = [checkedRxns;allRxnsWithMet];
    for i = allRxnsWithMet
        if ~ismember(i,checkedRxns) %don't check rxns we already checked
            %check all the mets in the reaction
            allMetsInRxn = find(model.S(:,i) ~=0);
            if length(allMetsInRxn) > 2
               b = false;
               return;
            end
            for j = allMetsInRxn
                if ~checkMetRec(model, checkedRxnsNextLevel, j)
                    b=false;
                    return;
                end
            end
        end
    end
end
%}
