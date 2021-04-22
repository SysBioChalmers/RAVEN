function [reducedModel, removedRxns, indexedDuplicateRxns]=contractModel(model,distReverse)
% contractModel
%   Contracts a model by grouping all identical reactions. Similar to the
%   deleteDuplicates part in simplifyModel but more care is taken here
%   when it comes to gene associations
%
%   model                  a model structure
%   distReverse            distinguish reactions with same metabolites
%                          but different reversibility as different
%                          reactions (opt, default true)
%
%   reducedModel           a model structure without duplicate reactions
%   removedRxns            cell array for the removed duplicate reactions
%   indexedDuplicateRxns   indexed cell array for the removed duplicate
%                          reactions (multiple valuess separated by semicolon)
%
%   NOTE: This code might not work for advanced grRules strings
%         that involve nested expressions of 'and' and 'or'.
%
%   Usage: [reducedModel, removedRxns, indexedDuplicateRxns]=contractModel(model,distReverse)

if nargin<2
    distReverse=true;
end

%First sort the model so that reversible reactions are in the same
%direction
modelS=sortModel(model);

%Get a list of duplicate reactions
if distReverse
    x=[modelS.S; model.rev']';
else
    x=modelS.S';
end
[~,I,J] = unique(x,'rows','first');

%Initialize cell array of indexedDuplicateRxns
indexedDuplicateRxns=cell(numel(model.rxns),1);
indexedDuplicateRxns(:)={''};

duplicateRxns=setdiff(1:numel(model.rxns),I);
mergeTo=I(J(duplicateRxns));

mergedRxns=unique(mergeTo);

%Now add all the info from this one. Print a warning if they have different
%bounds or objective function coefficients. Uses the widest bounds and
%largest magnitude of objective coefficient
for i=1:numel(mergedRxns)
    duplRxn=transpose([mergedRxns(i),duplicateRxns(mergeTo==mergedRxns(i))]);
    if numel(unique(model.lb(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' have different lower bound. Uses the most negative/smallest lower bound'];
        dispEM(EM,false);
        model.lb(mergedRxns(i))=min(model.lb(duplRxn));
    end
    if numel(unique(model.ub(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' have different upper bound. Uses the most positive/largest upper bound'];
        dispEM(EM,false);
        model.ub(mergedRxns(i))=max(model.ub(duplRxn));
    end
    if numel(unique(model.c(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' has a different objective function coefficient. Uses the largest coefficient'];
        dispEM(EM,false);
        model.c(mergedRxns(i))=max(model.c(duplRxn));
    end
    if isfield(model,'grRules') && any(~isempty(model.grRules(duplRxn)))
        rules=model.grRules(duplRxn);
        allRules={};
        for j=1:numel(rules)
            rules{j}=ignoreORinComplex(rules{j});
            allRules=[allRules regexp(rules{j},' or ','split')];
        end
        allRules=unique(allRules);
        allRules=strrep(allRules,'__OR__',' or ');
        andRules=~isempty(strfind(allRules,' and '));
        allRules(andRules)=strcat('(',allRules(andRules),')');
        if numel(allRules)==1
            model.grRules{mergedRxns(i)}=allRules{1};
        else
            model.grRules{mergedRxns(i)}=strjoin(allRules,' or ');
        end
    end    
    if isfield(model,'eccodes') && any(~isempty(model.eccodes(duplRxn)))
        codes=model.eccodes(duplRxn);
        allCodes={};
        for j=1:numel(codes)
            allCodes=[allCodes regexp(codes{j},';','split')];
        end
        allCodes=unique(allCodes);
        if numel(allCodes)==1
            model.eccodes{mergedRxns(i)}=allCodes{1};
        else
            model.eccodes{mergedRxns(i)}=strjoin(allCodes,';');
        end
    end
    %Generate indexedDuplicateRxns cell array
    if numel(duplRxn)==2
        indexedDuplicateRxns{duplRxn(1)}=model.rxns{duplRxn(2)};
    else
        indexedDuplicateRxns{duplRxn(1)}=strjoin(model.rxns(duplRxn(2:end)),';');
    end
end

%Delete the duplicate reactions
reducedModel=removeReactions(model,duplicateRxns);
removedRxns=model.rxns(duplicateRxns);
[~, index]=ismember(reducedModel.rxns,model.rxns);
indexedDuplicateRxns=indexedDuplicateRxns(index);

if isfield(reducedModel,'rxnGeneMat')
    %Fix grRules and reconstruct rxnGeneMat
    [grRules,rxnGeneMat] = standardizeGrRules(reducedModel,true);
    reducedModel.grRules = grRules;
    reducedModel.rxnGeneMat = rxnGeneMat;
end
end

function grRule = ignoreORinComplex(grRule)
%In a grRule, if OR relationship is nested in an AND relationship, then
%obfuscate the OR before splitting isoenzymes
grRule=['(' grRule ')'];
brOpen=strfind(grRule,'(');
brClose=strfind(grRule,')');
andPos=strfind(grRule,' and ');
%Find opening bracket closest before AND
stillCapturing = 0;
for i=1:numel(andPos)
    searchPos = andPos(i);
    while stillCapturing == 0
        closestOpen = brOpen(max(find(brOpen<searchPos)));
        inbetweenClose = brClose(brClose<searchPos & brClose>closestOpen);
        if ~isempty(inbetweenClose)
            searchPos=max(inbetweenClose);
        else
            stillCapturing = 1;
            beginPos = closestOpen;
        end
    end
    stillCapturing = 0;
    searchPos = andPos(i);
    while stillCapturing == 0
        closestClose = brClose(min(find(brClose>searchPos)));
        inbetweenOpen = brOpen(brOpen>searchPos & brOpen<closestOpen);
        if ~isempty(inbetweenOpen)
            searchPos=min(closestClose);
        else
            stillCapturing = 1;
            endPos = closestClose;
        end
    end
    replacePart=regexprep(grRule(beginPos:endPos),' or ','__OR__');
    grRule=strrep(grRule,grRule(beginPos:endPos),replacePart);
end
grRule=grRule(2:end-1);
end
