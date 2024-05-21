function model = changeGeneAssoc(model,rxnID,geneAssoc,replace)
% changeGeneAssoc
%   Changes gene association for a defined reaction
%
%   model       a model structure to change the gene association
%   rxnID       string or cell array of reaction IDs
%   geneAssoc   string of additional or replacement gene association.
%               Should be written with ' and ' to indicate subunits, ' or '
%               to indicate isoenzymes, and brackets '()' to separate
%               different instances
%   replace     true if old gene association should be replaced with new
%               association. False if new gene association should be
%               concatenated to the old association (optional, default true)
%
%   model       an updated model structure
%
% Usage: changeGeneAssoc(model,rxnID,geneAssoc,replace)

if nargin==3
    replace=true;
end

if isstr(rxnID)
    rxnID={rxnID};
end

if ~isstr(geneAssoc)
    EM='geneAssoc has to be string';
    error(sprintf(EM))
end


% Add genes to model
geneList=transpose(cell(unique(regexp(geneAssoc,'[)(]*|( and )*|( or )*','split')))); % Extract individual, unique genes from the geneAssoc provided
geneList=geneList(~cellfun(@isempty, geneList));
genesToAdd.genes=setdiff(geneList,model.genes); % Only keep the genes that are not yet part of the model.genes.
if ~isempty(genesToAdd.genes)
    model=addGenesRaven(model,genesToAdd); % Add genes
end

% Find reaction and gene indices
idx=getIndexes(model,rxnID,'rxns');

% Change gene associations
if replace==true % Replace old gene associations
    model.grRules(idx)=cellstr(geneAssoc);
else % Add gene associations, add new gene rules after 'OR'.
    model.grRules(idx)=cellstr('(',[model.grRules{idx},') or (',geneAssoc,')']);
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
end
