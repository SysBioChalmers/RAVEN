function model = changeGrRules(model,rxns,grRules,replace)
% changeGrRules
%   Changes multiple grRules at the same time.
%
%   model       a model structure to change the gene association
%   rxns        string or cell array of reaction IDs
%   grRules     string of additional or replacement gene association.
%               Should be written with ' and ' to indicate subunits, ' or '
%               to indicate isoenzymes, and brackets '()' to separate
%               different instances
%   replace     true if old gene association should be replaced with new
%               association. False if new gene association should be
%               concatenated to the old association (opt, default true)
%
%   model       an updated model structure
%
%   Usage: changeGrRules(model,rxns,grRules,replace)

if nargin==3
    replace=true;
end

rxns=convertCharArray(rxns);
grRules=convertCharArray(grRules);

if ~(numel(grRules)==numel(rxns))
    error('Number of rxns and grRules should be identical')
end

for i=1:length(rxns)
    % Add genes to model
    geneList=transpose(cell(unique(regexp(grRules{i},'[)(]*|( and )*|( or )*','split')))); % Extract individual, unique genes from the geneAssoc provided
    geneList=geneList(~cellfun(@isempty, geneList));
    genesToAdd.genes=setdiff(geneList,model.genes); % Only keep the genes that are not yet part of the model.genes.
    if ~isempty(genesToAdd.genes)
        model=addGenesRaven(model,genesToAdd); % Add genes
    end
    
    % Find reaction and gene indices
    idx=getIndexes(model,rxns,'rxns');
end

% Change gene associations
if replace==true % Replace old gene associations
    model.grRules(idx)=grRules;
else % Add gene associations, add new gene rules after 'OR'.
    model.grRules(idx)=strcat('(',model.grRules(idx),') or (',grRules,')');
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
end
