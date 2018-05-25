function rxnGeneMat=getRxnGeneMat(model)
% getRxnGeneMat
%   Re-generate rxnGeneMat matrix from genes and grRules fields
%
%   model          a RAVEN model structure
%
%   Usage: rxnGeneMat=getRxnGeneMat(model)
%
%   Simonas Marcisauskas, 2017-12-04
%


%Check gene association for each reaction and populate rxnGeneMat
if ~isempty(model.genes)
    rxnGeneMat=zeros(numel(model.rxns),numel(model.genes));
end
if ~isempty(model.grRules)
    tempRules=model.grRules;
    for i=1:length(model.rxns)
       %Check that all gene associations have a match in the gene list
       if ~isempty(model.grRules{i})          
           tempRules{i}=regexprep(tempRules{i},' and | or ','>'); %New format: Genes are separated 'and' and 'or' strings with parentheses
           tempRules{i}=regexprep(tempRules{i},'(',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
           tempRules{i}=regexprep(tempRules{i},')',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
           indexesNew=strfind(tempRules{i},'>'); %Old format: Genes are separated by ":" for AND and ";" for OR
           indexes=strfind(tempRules{i},':'); %Old format: Genes are separated by ":" for AND and ";" for OR
           indexes=unique([indexesNew indexes strfind(tempRules{i},';')]);
           if isempty(indexes)
               %See if you have a match
               I=find(strcmp(tempRules{i},model.genes));
               rxnGeneMat(i,I)=1;
           else
               temp=[0 indexes numel(tempRules{i})+1];
               for j=1:numel(indexes)+1
                   %The reaction has several associated genes
                   geneName=tempRules{i}(temp(j)+1:temp(j+1)-1);
                   I=find(strcmp(geneName,model.genes));
                   rxnGeneMat(i,I)=1;
               end
           end
       end
    end
end
rxnGeneMat=sparse(rxnGeneMat);

end
