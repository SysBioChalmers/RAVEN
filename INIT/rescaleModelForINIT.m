function model=rescaleModelForINIT(model, maxStoichDiff)
% rescaleModelForINIT
%
% The idea with this function is to rescale the MILP problem in ftINIT to avoid large differences
% in flux magnitudes between reactions. Such differences cause among other things
% difficulties regarding tolerances for integer variables.
% For now it just scales down all reactions with high stoichiometric coefficients
% There is room for improvement here - the best would be to convert mets such as albumin
% to instead represent 1/100 albumin - that would create much less extreme coefficients.
% This type of improvement is known as scaling in the literature around LPs and MILPs.
%
% model         the model to be modified (input and output)
% maxStoichVal  all reactions with stoichiometric coefficent higher than this 
%               will be scaled down. (opt, default 250)

if (nargin < 2)
    maxStoichDiff = 25;
end

%Define maxMinRatio to be the ratio between the largest and smallest 
%stoichiometric coefficients in each reaction.
%Scale all rxns with maxMinRatio > maxStoichDiff - just set all coeffs that are
%larger than maxStoichDiff*minVal to that value, and center the mean coeff to 1 for all rxns.
SAbs = abs(model.S);
for i = 1:numel(model.rxns)
    tmp = SAbs(:,i);
    tmp = tmp(tmp ~= 0);
    mn = min(tmp);
    %modify matrix
    sign = ones(numel(model.mets),1);
    sign(model.S(:,i) < 0) = -1;
    sel = SAbs(:,i) > maxStoichDiff*mn;
    model.S(sel,i) = sign(sel) .* maxStoichDiff*mn;
end
%center around 1
boolMat = model.S ~= 0;
absMat = abs(model.S);
rxnScales = sum(boolMat,1)./sum(absMat,1);

model.S = model.S .* rxnScales;

%{
if (nargin < 2)
    maxStoichVal = 250;
end


%find all reactions containing high stoichiometry coefficients
maxCoeffs = max(model.S, [], 1);
rxnsToCheck = maxCoeffs > maxStoichVal;

%for debugging:
%constructEquations(model, model.rxns(rxnsToCheck))

rxnsToCheckInd = find(rxnsToCheck);

%We just scale all reactions so they don't have a coefficient > maxStoichVal for now

%do this with a loop, doesn't really matter
for rxnInd = rxnsToCheckInd
    %find the lowest index which is not a pointless metabolite, such as H+, H2O etc.
    rxnS = model.S(:,rxnInd);
    %metInd = find(rxnS ~= 0);
    
    largestCoeff = max(abs(rxnS));
    scaleFactor = maxStoichVal/largestCoeff;
    model.S(:,rxnInd) = model.S(:,rxnInd) .* scaleFactor;
end
%}
end
