function balanceStructure=getElementalBalance(model,varargin)
% getElementalBalance  Check whether the reactions of a model are balanced.
%
% Checks a model to see if the reactions are elementally balanced.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% rxns : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of
%     indexes. Only these reactions will be checked (default model.rxns).
% printUnbalanced : logical
%     print warnings about the reactions that were unbalanced
%     (default false).
% printUnparsable : logical
%     print warnings about the reactions that cannot be parsed
%     (default false).
%
% Returns
% -------
% balanceStructure : struct
%     elemental balance structure with fields:
%
%     - balanceStatus : 1 if the reaction is balanced, 0 if it is
%       unbalanced, -1 if it could not be balanced due to missing
%       information, -2 if it could not be balanced due to an error.
%       Elemental only; charge is reported separately below
%     - chargeStatus : 1 if the reaction is charge balanced, 0 if it is
%       not, -1 if any participating metabolite has no charge (or the
%       model has no metCharges field), in which case the charge balance
%       is unknown rather than zero
%     - chargeResidual : the sum of charges over the reaction, NaN where
%       chargeStatus is -1
%     - elements : struct with fields abbrevs (cell array with
%       abbreviations for all used elements) and names (cell array with
%       the names for all used elements)
%     - leftComp : MxN matrix with the sum of coefficients for each of the
%       elements (N) for the left side of the reactions (M)
%     - rightComp : the corresponding matrix for the right side
%
% Examples
% --------
%     balanceStructure = getElementalBalance(model, rxns, printUnbalanced, printUnparsable);

p=parseRAVENargs(varargin, {'rxns',[]; 'printUnbalanced',false; 'printUnparsable',false});
rxns=p.rxns;
if ~isempty(rxns) && ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
end
printUnbalanced=p.printUnbalanced;
printUnparsable=p.printUnparsable;

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
        error('RAVEN:badInput', '%s', EM);
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
            %Do not change a more serious error to a less serious one
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
%done even if the reaction as a whole could not be balanced
total=abs(balanceStructure.rightComp-balanceStructure.leftComp)>10^-8; %To deal with roundoff error

%Get which reactions are unbalanced. Do not change an error to just
%unbalanced
balanceStructure.balanceStatus(any(total,2))=min(balanceStructure.balanceStatus(any(total,2)),0);

%Reactions with no metabolites are not balanced (empty S column skips both loops above)
emptyRxns = full(sum(model.S ~= 0, 1).' == 0);
balanceStructure.balanceStatus(emptyRxns) = min(-1, balanceStructure.balanceStatus(emptyRxns));

%The remaining ones are all balanced
balanceStructure.balanceStatus(isnan(balanceStructure.balanceStatus))=1;

%Charge balance. This is reported separately from balanceStatus, which
%callers such as removeBadRxns and printModelStats read as a purely
%elemental verdict.
balanceStructure.chargeStatus=zeros(numel(model.rxns),1);
balanceStructure.chargeResidual=nan(numel(model.rxns),1);
if ~isfield(model,'metCharges')
    balanceStructure.chargeStatus(:)=-1;
else
    for j=1:numel(model.rxns)
        %Only the participating metabolites may be touched: S is sparse and
        %0*NaN is NaN, so a single unset charge anywhere in the model would
        %otherwise poison every reaction. Summing with 'omitnan' instead
        %would be worse, silently reporting an unknown residual as 0.
        idx=find(model.S(:,j));
        if isempty(idx) || any(isnan(model.metCharges(idx)))
            balanceStructure.chargeStatus(j)=-1;
        else
            balanceStructure.chargeResidual(j)=sum(full(model.S(idx,j)).*model.metCharges(idx));
            if abs(balanceStructure.chargeResidual(j))>10^-8 %Roundoff error
                balanceStructure.chargeStatus(j)=0;
            else
                balanceStructure.chargeStatus(j)=1;
            end
        end
    end
end

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
            warning('RAVEN:warning', '%s', EM);
        else
            EM=['The reaction ' model.rxns{toPrint(i)} ' could not be balanced due to a parsing error'];
            warning('RAVEN:warning', '%s', EM);
        end
    else
        %Find the compounds that are not balanced
        notBalanced=find(total(toPrint(i),:));
        for j=1:numel(notBalanced)
            EM=['The reaction ' model.rxns{toPrint(i)} ' is not balanced with respect to ' balanceStructure.elements.names{notBalanced(j)}];
            warning('RAVEN:warning', '%s', EM);
        end
    end
end

% Re-order the structure entries so they are consistent with the ordering of
% the input reaction indexes
if ~isempty(rxns)
    rxns = getIndexes(model,rxns,'rxns');
    [~,i] = sort(rxns);
    balanceStructure.balanceStatus(i) = balanceStructure.balanceStatus;
    balanceStructure.chargeStatus(i) = balanceStructure.chargeStatus;
    balanceStructure.chargeResidual(i) = balanceStructure.chargeResidual;
    balanceStructure.leftComp(i,:) = balanceStructure.leftComp;
    balanceStructure.rightComp(i,:) = balanceStructure.rightComp;
end
