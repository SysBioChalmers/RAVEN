function balanceStructure=getElementalBalance(model,rxns,printUnbalanced,printUnparsable)
% getElementalBalance
%   Checks a model to see if the reactions are elementally balanced
%
%   model             a model structure
%   rxns              either a cell array of reaction IDs, a logical vector
%                     with the same number of elements as reactions in the model,
%                     of a vector of indexes. Only these reactions will be
%                     checked (opt, default model.rxns)
%   printUnbalanced   print warnings about the reactions that were
%                     unbalanced (opt, default false)
%   printUnparsable   print warnings about the reactions that cannot be
%                     parsed (opt, default false)
%
%   balanceStructure
%       balanceStatus   1 if the reaction is balanced, 0 if it's unbalanced,
%                      -1 if it couldn't be balanced due to missing information,
%                      -2 if it couldn't be balanced due to an error
%       elements
%           abbrevs     cell array with abbreviations for all used elements
%           names       cell array with the names for all used elements
%       leftComp        MxN matrix with the sum of coefficients for each of
%                       the elements (N) for the left side of the
%                       reactions (M)
%       rightComp       the corresponding matrix for the right side
%
%   Usage: balanceStructure=getElementalBalance(model,rxns,printUnbalanced,printUnparsable)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<2
    rxns=[];
end

if nargin<3
    printUnbalanced=false;
end

if nargin<4
    printUnparsable=false;
end

if ~isempty(rxns)
    indexes=~getIndexes(model,rxns,'rxns',true);
    model=removeReactions(model,indexes,true);
end

balanceStructure.balanceStatus=nan(numel(model.rxns),1);

%Get the formulas
if isfield(model,'metFormulas')
    [balanceStructure.elements, useMat, exitFlag]=parseFormulas(model.metFormulas, true);
else
    if isfield(model,'inchis')
        [balanceStructure.elements, useMat, exitFlag]=parseFormulas(model.inchis, true,true);
    else
        EM='The model must contain either the "metFormulas" or the "inchis" field in order to test for elemental balancing';
        dispEM(EM);
    end
end

balanceStructure.leftComp=zeros(numel(model.rxns),numel(balanceStructure.elements.names));
balanceStructure.rightComp=balanceStructure.leftComp;

%Look at the right and left side of the reactions separately
S{1}=model.S;
S{2}=model.S;
S{1}(S{1}>0)=0; %Left side
S{2}(S{2}<0)=0; %Right side
S{1}=abs(S{1}); %Both should have positive values

%First do it for left side and then for right
for i=1:2
    for j=1:numel(model.rxns)
        %Get the balance status of the involved mets
        I=exitFlag(S{i}(:,j)~=0);
        if any(I==-1)
            balanceStructure.balanceStatus(j)=-2;
        end
        if any(I==0)
            %Don't change a more serious error to a less serious one
            balanceStructure.balanceStatus(j)=min(-1,balanceStructure.balanceStatus(j));
        end
        %Loop through each element
        for k=1:numel(balanceStructure.elements.names)
            if i==1
                balanceStructure.leftComp(j,k)=sum(S{i}(:,j).*useMat(:,k));
            else
                balanceStructure.rightComp(j,k)=sum(S{i}(:,j).*useMat(:,k));
            end
        end
    end
end

%Now compare the left and right sides to find which are unbalanced. This is
%done even if the reaction as a whole couldn't be balanced
total=abs(balanceStructure.rightComp-balanceStructure.leftComp)>10^-8; %To deal with roundoff error

%Get which reactions are unbalanced. Don't change an error to just
%unbalanced
balanceStructure.balanceStatus(any(total,2))=min(balanceStructure.balanceStatus(any(total,2)),0);

%The remaining ones are all balanced
balanceStructure.balanceStatus(isnan(balanceStructure.balanceStatus))=1;

%Print warnings
toPrint=[];
if printUnbalanced==true
    toPrint=[toPrint;find(balanceStructure.balanceStatus==0)];
end
if printUnparsable==true
    toPrint=[toPrint;find(balanceStructure.balanceStatus<0)];
end

toPrint=sort(toPrint);
for i=1:numel(toPrint)
    if balanceStructure.balanceStatus(toPrint(i))<0
        if balanceStructure.balanceStatus(toPrint(i))==-1
            EM=['The reaction ' model.rxns{toPrint(i)} ' could not be balanced due to missing information'];
            dispEM(EM,false);
        else
            EM=['The reaction ' model.rxns{toPrint(i)} ' could not be balanced due to a parsing error'];
            dispEM(EM,false);
        end
    else
        %Find the compounds that it's not balanced for
        notBalanced=find(total(toPrint(i),:));
        for j=1:numel(notBalanced)
            EM=['The reaction ' model.rxns{toPrint(i)} ' is not balanced with respect to ' balanceStructure.elements.names{notBalanced(j)}];
            dispEM(EM,false);
        end
    end
end
end
