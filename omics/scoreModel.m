function [rxnScores, geneScores, hpaScores, arrayScores]=scoreModel(model,hpaData,varargin)
% scoreModel  Score the reactions and genes in a model, tINIT-style.
%
% Scores reactions and genes in a model based on the expression data from HPA
% and/or gene expression data. Only used by getINITModel; ftINIT calls
% scoreComplexModel directly.
%
% This is a wrapper over scoreComplexModel, keeping this function's argument
% order and its convention of scoring an unmeasured gene as -Inf. The scoring
% itself is scoreComplexModel's, with the reduction over each grRule fixed to
% the operator this function's multipleGeneScoring selects.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% hpaData : struct
%     HPA data structure from parseHPA (optional if arrayData is supplied,
%     default []).
%
% Name-Value Arguments
% --------------------
% arrayData : struct
%     gene expression data structure (optional if hpaData is supplied,
%     default []). With fields:
%
%     - genes : cell array with the unique gene names.
%     - tissues : cell array with the tissue names. The list may not be
%       unique, as there can be multiple cell types per tissue.
%     - celltypes : cell array with the cell type names for each tissue.
%     - levels : GENESxTISSUES array with the expression level for each gene
%       in each tissue/celltype. NaN should be used when no measurement was
%       performed.
%     - threshold : a single value or a vector of gene expression thresholds,
%       above which genes are considered to be "expressed" (optional).
% tissue : char
%     tissue to score for. Should exist in either hpaData.tissues or
%     arrayData.tissues.
% celltype : char
%     cell type to score for. Should exist in either hpaData.celltypes or
%     arrayData.celltypes for this tissue (default is to use the maximum
%     values among all the cell types for the tissue).
% noGeneScore : double
%     score for reactions without genes (default -2).
% multipleGeneScoring : char
%     how to score reactions with several genes, "best" or "average"
%     (default "best"). "best" takes the highest-scoring gene, and is what
%     every caller uses. "average" averages down the grRule rather than flatly
%     over the reaction's genes, so a complex counts once against its
%     isozymes rather than once per subunit.
% multipleCellScoring : char
%     how to score when several cell types are used, "best" or "average"
%     (default "best").
% hpaLevelScores : struct
%     structure with numerical scores for the expression level categories
%     from HPA. The structure should have a "names" and a "scores" field
%     (default, see scoreComplexModel for the default scores).
%
% Returns
% -------
% rxnScores : double
%     scores for each of the reactions in model.
% geneScores : double
%     scores for each of the genes in model. Genes which are not in the
%     dataset(s) have -Inf as scores.
% hpaScores : double
%     scores for each of the genes in model if only taking hpaData into
%     account. Genes which are not in the dataset(s) have -Inf as scores.
% arrayScores : double
%     scores for each of the genes in model if only taking arrayData into
%     account. Genes which are not in the dataset(s) have -Inf as scores.
%
% Examples
% --------
%     [rxnScores, geneScores, hpaScores, arrayScores] = scoreModel(model, ...
%         hpaData, arrayData, tissue, celltype, noGeneScore, ...
%         multipleGeneScoring, multipleCellScoring, hpaLevelScores);
%
% See Also
% --------
% scoreComplexModel, getINITModel

p=parseRAVENargs(varargin, {'arrayData',[]; 'tissue',[]; 'celltype',[]; 'noGeneScore',-2; 'multipleGeneScoring','best'; 'multipleCellScoring','best'; 'hpaLevelScores',[]});
arrayData=p.arrayData;
tissue=p.tissue;
celltype=p.celltype;
noGeneScore=p.noGeneScore;
multipleGeneScoring=p.multipleGeneScoring;
multipleCellScoring=p.multipleCellScoring;
hpaLevelScores=p.hpaLevelScores;

if isempty(multipleGeneScoring)
    multipleGeneScoring='best';
end
multipleGeneScoring=lower(char(multipleGeneScoring));
if ~ismember(multipleGeneScoring,{'best','average'})
    EM='Valid options for multipleGeneScoring are "best" or "average"';
    error('RAVEN:badInput', '%s', EM);
end
if isempty(multipleCellScoring)
    multipleCellScoring='best';
end
multipleCellScoring=lower(char(multipleCellScoring));
if ~ismember(multipleCellScoring,{'best','average'})
    EM='Valid options for multipleCellScoring are "best" or "average"';
    error('RAVEN:badInput', '%s', EM);
end

%One operator for both and/or, so a rule reduces to a single statistic over
%its genes, which is what scoring by rxnGeneMat amounts to
if strcmp(multipleGeneScoring,'best')
    ruleScoring='max';
else
    ruleScoring='average';
end
%This function calls the reduction over cell types "best" where
%scoreComplexModel calls it "max"
if strcmp(multipleCellScoring,'best')
    cellScoring='max';
else
    cellScoring='average';
end

[rxnScores, geneScores, hpaScores, arrayScores] = scoreComplexModel(model, ...
    hpaData, arrayData, tissue, ...
    'celltype', celltype, ...
    'noGeneScore', noGeneScore, ...
    'isozymeScoring', ruleScoring, ...
    'complexScoring', ruleScoring, ...
    'multipleCellScoring', cellScoring, ...
    'hpaLevelScores', hpaLevelScores, ...
    'dataPrecedence', 'reaction');

%scoreComplexModel marks an unmeasured gene NaN, because removeLowScoreGenes
%treats NaN as "no evidence" and leaves such genes in place. This function's
%callers test with isinf instead
geneScores(isnan(geneScores)) = -Inf;
end
