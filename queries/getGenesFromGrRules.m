function [genes,rxnGeneMat] = getGenesFromGrRules(grRules, varargin)
% getGenesFromGrRules  Extract gene list and rxnGeneMat from grRules array.
%
% Parameters
% ----------
% grRules : cell
%     a cell array of model grRules, from which a list of genes is to be
%     extracted. NOTE: Boolean operators can be text ("and", "or") or
%     symbolic ("&", "|"), but there must be a space between operators and
%     gene names/IDs.
%
% Name-Value Arguments
% --------------------
% originalGenes : cell
%     the original gene list from the model as reference.
%
% Returns
% -------
% genes : cell
%     a unique list of all gene IDs that appear in grRules.
% rxnGeneMat : double
%     (optional) a binary matrix indicating which genes participate in each
%     reaction, where rows correspond to reactions (entries in grRules) and
%     columns correspond to genes.
%
% Examples
% --------
%     [genes, rxnGeneMat] = getGenesFromGrRules(grRules, originalGenes);


% handle input arguments
p=parseRAVENargs(varargin, {'originalGenes',[]});
originalGenes=p.originalGenes;

% check if the grRules use written or symbolic boolean operators
if any(contains(grRules,{' & ',' | '}))
    % fix some potential missing spaces between parentheses and &/|
    grRules = regexprep(grRules,'\)&',') &');   % ")&"  ->  ") &"
    grRules = regexprep(grRules,'&\(','& (');   % "&("  ->  "& ("
    grRules = regexprep(grRules,'\)\|',') |');  % ")|"  ->  ") |"
    grRules = regexprep(grRules,'\|\(','| (');  % "|("  ->  "| ("
else
    % fix some potential missing spaces between parentheses and AND/OR
    grRules = regexprep(grRules,'\)and',') and');  % ")and" ->  ") and"
    grRules = regexprep(grRules,'and\(','and (');  % "and(" ->  "and ("
    grRules = regexprep(grRules,'\)or',') or');    % ")or"  ->  ") or"
    grRules = regexprep(grRules,'or\(','or (');    % "or("  ->  "or ("
    
    % convert "and" to "&" and "or" to "|" (easier to work with symbols)
    grRules = regexprep(grRules, ' or ', ' | ');
    grRules = regexprep(grRules, ' and ', ' & ');
end

% extract list of genes from each reaction
rxnGenes = cellfun(@(r) regexprep(unique(strsplit(r,{' | ',' & '})),'[\(\) ]+',''),grRules,'UniformOutput',false);

% construct new gene list
nonEmpty = ~cellfun(@isempty,rxnGenes);
genes = unique(transpose([rxnGenes{nonEmpty}]));
genes(cellfun(@isempty,genes)) = [];

if ~isempty(originalGenes)
    if ~isequal(sort(originalGenes), sort(genes))
        error('The grRules and original gene list are inconsistent!');
    else
        genes = originalGenes;
    end
end

% construct new rxnGeneMat (if requested)
if nargout > 1
    rxnGeneCell = cellfun(@(rg) ismember(genes,rg),rxnGenes,'UniformOutput',false);
    rxnGeneMat = sparse(double(horzcat(rxnGeneCell{:})'));
end