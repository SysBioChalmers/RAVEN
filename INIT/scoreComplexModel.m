function [rxnScores, geneScores, hpaScores, arrayScores] = scoreComplexModel(model,hpaData,arrayData,tissue,varargin)
% scoreComplexModel  Score reactions and genes in a model with complex gene rules.
%
% Scores the reactions and genes in a model containing complex gene rules
% based on expression data from HPA and/or gene arrays.
%
% It is highly recommended that the model grRules are "cleaned" and verified by
% the "cleanModelGeneRules" function prior to scoring the model with
% scoreComplexModel.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% hpaData : struct
%     HPA data structure from parseHPA (optional if arrayData is supplied,
%     default []).
% arrayData : struct
%     gene expression data structure (optional if hpaData is supplied, default
%     []). With fields:
%
%     - genes : cell array with the unique gene names.
%     - tissues : cell array with the tissue names. The list may not be unique,
%       as there can be multiple cell types per tissue.
%     - celltypes : cell array with the cell type names for each tissue.
%     - levels : GENESxTISSUES array with the expression level for each gene in
%       each tissue/celltype. NaN should be used when no measurement was
%       performed.
%     - threshold : a single value or a vector of gene expression thresholds,
%       above which genes are considered to be "expressed". (optional, by
%       default, the mean expression levels of each gene across all tissues in
%       arrayData will be used as the threshold values).
% tissue : char
%     tissue to score for. Should exist in either hpaData.tissues or
%     arrayData.tissues.
%
% Name-Value Arguments
% --------------------
% celltype : char
%     cell type to score for. Should exist in either hpaData.celltypes or
%     arrayData.celltypes for this tissue (default is to use the max values
%     among all the cell types for the tissue).
% noGeneScore : double
%     score for reactions without genes (default -2).
% isozymeScoring : char
%     determines how scores are calculated for reactions with multiple genes
%     joined by "OR" expression(s) ("min", "max", "median", "average")
%     (default "max").
% complexScoring : char
%     determines how scores are calculated for reactions with multiple genes
%     joined by "AND" expression(s) ("min", "max", "median", "average")
%     (default "min").
% multipleCellScoring : char
%     determines how scores are calculated when several cell types are used
%     ("max" or "average") (default "max").
% hpaLevelScores : struct
%     structure with numerical scores for the expression level categories from
%     HPA. The structure should have a "names" and a "scores" field (default,
%     see code for default scores).
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


p=parseRAVENargs(varargin, {'celltype',[]; 'noGeneScore',[]; 'isozymeScoring',[]; 'complexScoring',[]; 'multipleCellScoring',[]; 'hpaLevelScores',[]});
celltype=p.celltype;
noGeneScore=p.noGeneScore;
isozymeScoring=p.isozymeScoring;
complexScoring=p.complexScoring;
multipleCellScoring=p.multipleCellScoring;
hpaLevelScores=p.hpaLevelScores;
if isempty(noGeneScore)
    noGeneScore = -2;
end
if isempty(isozymeScoring)
    isozymeScoring = 'max';
end
if isempty(complexScoring)
    complexScoring = 'min';
end
if isempty(multipleCellScoring)
    multipleCellScoring = 'max';
end
if isempty(hpaLevelScores)
    % The first four are for APE, the other ones for staining
    hpaLevelScores.names = {'High' 'Medium' 'Low' 'None' 'Strong' 'Moderate' 'Weak' 'Negative' 'Not detected'};
    hpaLevelScores.scores = [20 15 10 -8 20 15 10 -8 -8];
end

if isempty(hpaData) && isempty(arrayData)
    EM = 'Must supply hpaData, arrayData or both';
    error('RAVEN:badInput', '%s', EM);
end
if ~ismember(lower(isozymeScoring),{'min','max','median','average'})
    EM = 'Valid options for isozymeScoring are "min", "max", "median", or "average"';
    error('RAVEN:badInput', '%s', EM);
end
if ~ismember(lower(complexScoring),{'min','max','median','average'})
    EM = 'Valid options for complexScoring are "min", "max", "median", or "average"';
    error('RAVEN:badInput', '%s', EM);
end
if ~ismember(lower(multipleCellScoring),{'max','average'})
    EM = 'Valid options for multipleCellScoring are "max" or "average"';
    error('RAVEN:badInput', '%s', EM);
end


% Throw an error if array data for only one tissue is supplied without
% specifying threshold values
if ~isempty(arrayData)
    if numel(unique(arrayData.tissues)) < 2
        if ~isfield(arrayData,'threshold') || isempty(arrayData.threshold)
            EM = 'arrayData must contain measurements for at least two celltypes/tissues since the score is calculated based on the expression level compared to the overall average';
            error('RAVEN:badInput', '%s', EM);
        end
    end
end

% Process arrayData.threshold if necessary
if isfield(arrayData,'threshold') && (numel(arrayData.threshold) == 1)
    % if only a single gene threshold value is provided, then just
    % duplicate this value for all genes.
    arrayData.threshold = arrayData.threshold*ones(size(arrayData.genes));
end

% This is so the code can ignore which combination of input data is used
if isempty(arrayData)
    arrayData.genes={};
    arrayData.tissues={};
    arrayData.celltypes={};
    arrayData.levels=[];
    arrayData.threshold=[];
end
if isempty(hpaData)
    hpaData.genes={};
    hpaData.tissues={};
    hpaData.celltypes={};
    hpaData.levels={};
    hpaData.types={};
    hpaData.reliabilities={};
    hpaData.gene2Level=[];
    hpaData.gene2Type=[];
    hpaData.gene2Reliability=[];
end

% Check that the tissue exists
if ~ismember(upper(tissue),upper(hpaData.tissues)) && ~ismember(upper(tissue),upper(arrayData.tissues))
    EM = 'The tissue name does not match';
    error('RAVEN:badInput', '%s', EM);
end
if any(celltype)
    % Check that both data types has cell type defined if that is to be used
    if ~isfield(hpaData,'celltypes') || ~isfield(arrayData,'celltypes')
        EM = 'Both hpaData and arrayData must contain cell type information if cell type is to be used';
        error('RAVEN:badInput', '%s', EM);
    end
    if ~ismember(upper(celltype),upper(hpaData.celltypes)) && ~ismember(upper(celltype),upper(arrayData.celltypes))
        EM = 'The cell type name does not match';
        error('RAVEN:badInput', '%s', EM);
    end
end

% Some preprocessing of the structures to increase efficiency
% Remove all tissues that are not the correct one
J = ~strcmpi(hpaData.tissues,tissue);

%If cell type is supplied, then only keep that cell type
if any(celltype)
    J = J | ~strcmpi(hpaData.celltypes,celltype);
end

hpaData.tissues(J)=[];
if isfield(hpaData,'celltypes')
    hpaData.celltypes(J)=[];
end
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(:,J)=[];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(:,J)=[];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(:,J)=[];
end

% Remove all genes from the structures that are not in model or that are not
% measured in the tissue
if ~isempty(hpaData.genes) % This should not be necessary, but the summation is a 0x1 matrix and the other is []
    I = ~ismember(hpaData.genes,model.genes) | sum(hpaData.gene2Level,2) == 0;
else
    I = [];
end
hpaData.genes(I) = [];
if isfield(hpaData,'gene2Level')
    hpaData.gene2Level(I,:) = [];
end
if isfield(hpaData,'gene2Type')
    hpaData.gene2Type(I,:) = [];
end
if isfield(hpaData,'gene2Reliability')
    hpaData.gene2Reliability(I,:) = [];
end

I = strcmpi(arrayData.tissues,tissue);
% If cell type is supplied, then only keep that cell type
if any(celltype)
    I = I & strcmpi(arrayData.celltypes,celltype);
end

% Remove all genes from the structures that are not in model or that are not
% measured in the tissue
J = ~ismember(arrayData.genes,model.genes) | all(isnan(arrayData.levels(:,I)),2);
arrayData.genes(J) = [];
arrayData.levels(J,:) = [];
if isfield(arrayData,'threshold')
    arrayData.threshold(J) = [];
end

% Calculate the scores for the arrayData. These scores are calculated for
% each genes from its fold change between the tissue/celltype(s) in question
% and all other celltypes, or the threshold if supplied. This is a lower 
% quality data than protein abundance, since gene abundance is an indirect 
% estimate of protein level. These scores are therefore only used for genes
% for which there is no HPA data available. The fold changes are transformed
% as min(5*log(x),10) for x > 1 and max(5*log(x),-5) for x < 1 in order to 
% have negative scores for lower expressed genes and to scale the scores
% to have somewhat lower weights than the HPA scores
tempArrayLevels = arrayData.levels;
tempArrayLevels(isnan(tempArrayLevels)) = 0;
if isfield(arrayData,'threshold') && ~isempty(arrayData.threshold)
    % if provided, the user-supplied expression threshold value(s) will be
    % used as the "average" expression level to which each gene is compared
    average = arrayData.threshold;
else
    average = sum(tempArrayLevels,2)./sum(~isnan(arrayData.levels),2);
end
if strcmpi(multipleCellScoring,'max')
    current = max(tempArrayLevels(:,I),[],2);
else
    current = sum(tempArrayLevels(:,I),2)./sum(~isnan(arrayData.levels(:,I)),2);
end
if ~isempty(current)
    aScores = 5*log(current./average);
else
    aScores = [];
end
aScores(aScores > 0) = min(aScores(aScores > 0),10);
aScores(aScores < 0) = max(aScores(aScores < 0),-5);
aScores(isnan(aScores)) = -5;  % NaNs occur when gene expression is zero across all tissues

% Map the HPA levels to scores
[I, J] = ismember(upper(hpaData.levels),upper(hpaLevelScores.names));
if ~all(I)
    EM = 'There are expression level categories that do not match to hpaLevelScores';
    error('RAVEN:badInput', '%s', EM);
end
[K, ~, M] = find(hpaData.gene2Level);
scores = hpaLevelScores.scores(J);
% Reduce over the cell types each gene was actually measured in. Building a
% sparse GENESxCELLTYPES matrix instead would put a numeric 0 at every
% unmeasured (gene, cell type) pair, and 0 sits between 'Low' (10) and
% 'Not detected' (-8): a gene not detected in one of two cell types would
% score 0 rather than -8 under 'max', and would never be pruned. The
% arrayData branch above likewise divides by the number of measurements.
measured = reshape(scores(M), [], 1);
if strcmpi(multipleCellScoring,'max')
    hScores = accumarray(K(:), measured, [numel(hpaData.genes) 1], @max, 0);
else
    hScores = accumarray(K(:), measured, [numel(hpaData.genes) 1], @mean, 0);
end

% Assign gene scores, prioritizing HPA (protein) data over arrayData (RNA)
geneScores = nan(numel(model.genes),1);
hpaScores = -Inf(numel(model.genes),1);
arrayScores = -Inf(numel(model.genes),1);

[I, J] = ismember(model.genes,arrayData.genes);
geneScores(I) = aScores(J(I));
arrayScores(I) = aScores(J(I));

[I, J] = ismember(model.genes,hpaData.genes);
geneScores(I) = hScores(J(I));
hpaScores(I) = hScores(J(I));


% To speed things up, only need to score each unique grRule once
[uRules,~,rule_ind] = unique(model.grRules);
uScores = nan(size(uRules));

% convert logic operators to symbols
uRules = regexprep(uRules,' and ',' & ');
uRules = regexprep(uRules,' or ',' | ');

% score based on presence/combination of & and | operators
for i = 1:numel(uRules)
    if isempty(uRules{i})
        uScores(i) = noGeneScore;
    elseif contains(uRules{i},'&') && contains(uRules{i},'|')
        uScores(i) = scoreComplexRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    else
        uScores(i) = scoreSimpleRule(uRules{i},model.genes,geneScores,isozymeScoring,complexScoring);
    end
end

% NaN reaction scores should be changed to the no-gene score
uScores(isnan(uScores)) = noGeneScore;

% re-map unique rule scores to model
rxnScores = uScores(rule_ind);

end



function rScore = scoreSimpleRule(rule,genes,gScores,isozymeScoring,complexScoring)
%Score reactions with simple gene rules (those with all ANDs or all ORs).

if ~contains(rule,'&')
    scoreMethod = isozymeScoring;
elseif ~contains(rule,'|')
    scoreMethod = complexScoring;
else
    error('This function cannot handle complex gene rules.');
end
ruleGenes = unique(regexp(rule,'[^&|\(\) ]+','match'));
[~,geneInd] = ismember(ruleGenes,genes);
switch lower(scoreMethod)
    case 'min'
        rScore = min(gScores(geneInd),[],'omitnan');
    case 'max'
        rScore = max(gScores(geneInd),[],'omitnan');
    case 'median'
        rScore = median(gScores(geneInd),'omitnan');
    case 'average'
        rScore = mean(gScores(geneInd),'omitnan');
end

end



function rScore = scoreComplexRule(rule,genes,gScores,isozymeScoring,complexScoring)
%Score reactions with complex gene rules (those with both ANDs or ORs).

% Specify phrases to search for in the grRule. These phrases will find
% genes grouped by all ANDs (first phrase) or all ORs (second phrase).
search_phrases = {'\([^&|\(\) ]+( & [^&|\(\) ]+)+\)', '\([^&|\(\) ]+( \| [^&|\(\) ]+)+\)'};

% initialize some variables
subsets = {};  % subsets are groups of genes grouped by all ANDs or all ORs
c = 1;  % counter to keep track of the group (subset) number
r_orig = rule;  % record original rule to determine when it stops changing
for k = 1:100  % iterate some arbitrarily high number of times
    for j = 1:length(search_phrases)
        new_subset = regexp(rule,search_phrases{j},'match')';  % extract subsets
        if ~isempty(new_subset)
            subsets = [subsets; new_subset];  % append to list of subsets
            subset_nums = arrayfun(@num2str,(c:length(subsets))','UniformOutput',false);  % get group numbers to be assigned to the new subsets, and convert to strings
            rule = regexprep(rule,search_phrases{j},strcat('#',subset_nums,'#'),'once');  % replace the subsets in the expression with their group numbers (enclosed by "#"s)
            c = c + length(new_subset);
        end
    end
    if isequal(rule,r_orig)
        break;  % stop iterating when rule stops changing
    else
        r_orig = rule;
    end
end
subsets{end+1} = rule;  % add final state of rule as the last subset

% score each subset and append to gene list and gene scores
for i = 1:numel(subsets)
    gScores = [gScores; scoreSimpleRule(subsets{i},genes,gScores,isozymeScoring,complexScoring)];
    genes = [genes; {strcat('#',num2str(i),'#')}];
end

% the final subset score is the overall reaction score
rScore = gScores(end);

end

