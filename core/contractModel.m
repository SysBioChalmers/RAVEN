function [reducedModel removedRxns]=contractModel(model)
% contractModel
%   Contracts a model by grouping all identical reactions. Similar to the
%   deleteDuplicates part in simplifyModel but more care is taken here
%   when it comes to gene associations
%
%	model           a model structure
%
%	reducedModel	a model structure with grouped reactions
%   removedRxns     cell array with the duplicate reactions
%
%   NOTE: This code might not work for advanced grRules strings
%         that involve nested expressions of 'and' and 'or'.
%
%   Usage: [reducedModel removedRxns]=contractModel(model)
%
%   Rasmus Agren, 2013-08-01
%

%First sort the model so that reversible reactions are in the same
%direction
modelS=sortModel(model);

%Get a list of duplicate reactions
x=[modelS.S; model.rev']';
[B,I,J] = unique(x,'rows','first');

duplicateRxns=setdiff(1:numel(model.rxns),I);
mergeTo=I(J(duplicateRxns));

%Now add all the info from this one. Print a warning if they have different
%bounds or objective function coefficients. Uses the widest bounds and largest
%magnitude of objective coefficient
for i=1:numel(duplicateRxns)
    if model.lb(duplicateRxns(i))<model.lb(mergeTo(i))
       dispEM(['Duplicate reaction ' model.rxns{duplicateRxns(i)} ' has wider lower bound. Uses the most negative/smallest lower bound'],false); 
       model.lb(mergeTo(i))=model.lb(duplicateRxns(i));
    end
    if model.ub(duplicateRxns(i))>model.ub(mergeTo(i))
       dispEM(['Duplicate reaction ' model.rxns{duplicateRxns(i)} ' has wider upper bound. Uses the most positive/largest upper bound'],false); 
       model.ub(mergeTo(i))=model.ub(duplicateRxns(i));
    end
    if abs(model.c(duplicateRxns(i)))>abs(model.c(mergeTo(i)))
       dispEM(['Duplicate reaction ' model.rxns{duplicateRxns(i)} ' has a larger objective function coefficient. Uses the largest coefficient'],false); 
       model.c(mergeTo(i))=model.c(duplicateRxns(i));
    end
    
    %Genes are added as 'or'
    if isfield(model,'rxnGeneMat')
        commonGenes=find(model.rxnGeneMat(duplicateRxns(i),:) & model.rxnGeneMat(mergeTo(i),:));
        newGenes=model.rxnGeneMat(duplicateRxns(i),:);
        newGenes(commonGenes)=0;
        model.rxnGeneMat(mergeTo(i),:)=model.rxnGeneMat(mergeTo(i),:)+newGenes;
    end
    if isfield(model,'grRules')
        if any(model.grRules{duplicateRxns(i)})
           if any(model.grRules{mergeTo(i)})
               %Split both grStrings on ' or ' and then put together union
               %with ' or '
               rules1=regexp(model.grRules{mergeTo(i)},' or ','split');
               rules2=regexp(model.grRules{duplicateRxns(i)},' or ','split');
               allRules=union(rules1,rules2);
               
               %Probably not the nicest way to do this
               model.grRules{mergeTo(i)}=allRules{1};
               for j=2:numel(allRules)
                   model.grRules{mergeTo(i)}=[model.grRules{mergeTo(i)} ' or ' allRules{j}];
               end
           else
               model.grRules{mergeTo(i)}=model.grRules{duplicateRxns(i)};
           end
        end
    end
    if isfield(model,'eccodes')
        if any(model.eccodes{duplicateRxns(i)})
           if any(model.eccodes{mergeTo(i)})
               %Split on ';' and put together the union with ';'
               codes1=regexp(model.eccodes{mergeTo(i)},';','split');
               codes2=regexp(model.eccodes{duplicateRxns(i)},';','split');
               codes=union(codes1,codes2);
               model.eccodes{mergeTo(i)}=codes{1};
               for j=2:numel(codes)
                  model.eccodes{mergeTo(i)}=[model.eccodes{mergeTo(i)} ';' codes{j}]; 
               end
           else
               model.eccodes{mergeTo(i)}=model.eccodes{duplicateRxns(i)};
           end
        end
    end
end

%Delete the duplicate reactions
reducedModel=removeReactions(model,duplicateRxns);
removedRxns=model.rxns(duplicateRxns);
end
