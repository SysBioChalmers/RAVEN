function model = changeGrRules(model,rxns,grRules,replace)
% changeGrRules  Change multiple grRules at the same time.
%
% Parameters
% ----------
% model : struct
%     a model structure to change the gene association.
% rxns : char or cell
%     reaction IDs.
% grRules : char or cell
%     additional or replacement gene association. Should be written with
%     ' and ' to indicate subunits, ' or ' to indicate isoenzymes, and
%     brackets '()' to separate different instances.
% replace : logical, optional
%     true if old gene association should be replaced with new association.
%     False if new gene association should be concatenated to the old
%     association (default true).
%
% Returns
% -------
% model : struct
%     an updated model structure.
%
% Examples
% --------
%     model = changeGrRules(model, rxns, grRules, replace);

if nargin==3
    replace=true;
end

rxns=convertCharArray(rxns);
grRules=convertCharArray(grRules);

if isscalar(grRules) && ~isscalar(rxns)
    grRules = repmat(grRules,1,numel(rxns));
end
if ~(numel(grRules)==numel(rxns))
    error('Number of rxns and grRules should be identical')
end

% Add genes to model
geneList = getGenesFromGrRules(grRules);
genesToAdd.genes=setdiff(geneList,model.genes); % Only keep the genes that are not yet part of the model.genes.
if ~isempty(genesToAdd.genes)
    model=addGenesRaven(model,genesToAdd); % Add genes
end
    
% Find reaction and gene indices
idx=getIndexes(model,rxns,'rxns');

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
