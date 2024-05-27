function expr = getExprForRxnScore(scores, threshold)
% getExprForRxnScore
%   Converts a reaction score to the gene expression (CPM or TPM) required 
%   to get that reaction score, if the GPR is only a single gene.
%   Useful function primarily in test cases, where you want to be able to
%   define the reaction scores of rxns, but need to send in gene expression.
%
%   scores        Vector of scores to convert
%   threshold     Gene threshold (optional, default 1)
%
%   expr          The resulting gene expression vector.
%
% Usage: expr = getExprForRxnScore(scores, threshold)

if nargin < 2
    threshold = 1;
end

%This is how the score is calculated: 5*log(expression./threshold)
%expression = threshold*10.^(scores/5)
%This is a bit confusing - it seems that it is threshold*e.^(scores/5)
%This is probably what is being used in scoreComplexModel, this code
%negates that perfectly, see the T0009 test case.

expr = threshold*exp(scores/5);

end
