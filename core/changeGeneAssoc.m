function model = changeGeneAssoc(model,rxnID,geneAssoc,replace)
% changeGeneAssoc
%   Changes gene association for a defined reaction
%
%   model       a model structure to change the gene association
%   rxnID       string of reaction ID
%   geneAssoc   string of additional or replacement gene association.
%               Should be written with ' and ' to indicate subunits, ' or '
%               to indicate isoenzymes, and brackets '()' to separate
%               different instances
%   replace     true if old gene association should be replaced with new
%               association. False if new gene association should be
%               concatenated to the old association (opt, default true)
%
%   model       an updated model structure
%
%   Usage: changeGeneAssoc(model,rxnID,geneAssoc,replace)
%
%   Eduard Kerkhoven, 2015-10-07

% 
if nargin==3
    replace=true;
end

% Add genes to model
geneList=transpose(cell(unique(regexp(geneAssoc,'[)(]*|( and )*|( or )*','split')))); % Extract individual, unique genes from the geneAssoc provided
genesToAdd.genes=setdiff(geneList(2:length(geneList)),model.genes); % Only keep the genes that are not yet part of the model.genes.
model=addGenes(model,genesToAdd); % Add genes

% Find reaction and gene indices
[~,idxRxn]=ismember(rxnID,model.rxns); % Find indices of rxnID in model.rxns
idxGene=zeros(length(geneList)-1,1); % Prepare an empty list to be filed with indices, and subsequently fill it.
for n=2:length(geneList)
    idxGene(n-1)=find(strcmp(model.genes,geneList(n)));
end

% Change gene associations
if replace==true % Replace old gene associations
    model.grRules(idxRxn)=cellstr(geneAssoc);
    model.rxnGeneMat(idxRxn,:)=0;
    model.rxnGeneMat(idxRxn,idxGene)=1;
else % Add gene associations, add new gene rules after 'OR'.
    model.grRules(idxRxn)=cellstr([model.grRules{idxRxn},' or ',geneAssoc]);
    model.rxnGeneMat(idxRxn,idxGene)=1;
end
end