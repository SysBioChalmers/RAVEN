function model = renameModelGenes(model, geneTable, fromCol, toCol)
% renameModelGenes  Replace gene identifiers in a RAVEN model.
%
% Updates model.genes and model.grRules by substituting the identifiers
% currently used in the model (fromCol) with new identifiers (toCol) from
% a mapping table produced by getGeneData or supplied manually. After
% renaming, model.rxnGeneMat is rebuilt automatically.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct. Must contain model.genes and model.grRules;
%     model.rxnGeneMat is rebuilt automatically.
% geneTable : table | char
%     Gene-mapping table. Can be supplied as:
%       (a) A MATLAB table variable already loaded in the workspace, or
%       (b) A path to a tab-delimited .tsv file — loaded automatically.
% fromCol : char
%     Column in geneTable whose values match the identifiers currently in
%     model.genes (e.g. 'locus_tag').
% toCol : char
%     Column in geneTable whose values will replace the current
%     identifiers (e.g. 'gene_name').
%
% Returns
% -------
% model : struct
%     Updated model struct with:
%       model.genes      — renamed gene list.
%       model.grRules    — gene-reaction rules using new identifiers.
%       model.rxnGeneMat — rebuilt sparse matrix matching the new gene list.
%
% Examples
% --------
%     model     = readYAMLmodel('iCre1355.yml');
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv');
%     model     = renameModelGenes(model, geneTable, 'locus_tag', 'gene_name');
%
%     % Supply a TSV file path directly
%     model = renameModelGenes(model, 'chlamy.tsv', 'locus_tag', 'gene_name');
%
% Notes
% -----
% * Genes with no entry in fromCol, or whose toCol value is empty, are left
%   unchanged; a warning lists them for investigation.
% * Word-boundary matching prevents partial substitution (e.g. 'gene1'
%   is not replaced inside 'gene10').
% * grRules are standardized via standardizeGrRules before and after
%   renaming to ensure consistent formatting.

    if nargin < 1 || isempty(model)
        error('renameModelGenes:missingInput', ...
            'model is required. Provide a RAVEN model struct.');
    end

    if nargin < 2 || isempty(geneTable)
        error('renameModelGenes:missingInput', ...
            'geneTable is required. Provide a table variable or a path to a .tsv file.');
    end

    if nargin < 3 || isempty(fromCol)
        error('renameModelGenes:missingInput', ...
            'fromCol is required. Specify the column in geneTable whose values match model.genes.');
    end

    if nargin < 4 || isempty(toCol)
        error('renameModelGenes:missingInput', ...
            'toCol is required. Specify the column in geneTable with the replacement identifiers.');
    end

    fromCol = char(fromCol);
    toCol   = char(toCol);

    % Load table if a file path was supplied
    if ischar(geneTable) || isstring(geneTable)
        if ~isfile(geneTable)
            error('renameModelGenes:fileNotFound', ...
                  'Table file not found: %s', geneTable);
        end
        geneTable = readtable(geneTable, 'Delimiter', '\t', 'FileType', 'text', ...
                              'TextType', 'char');
        fprintf('Loaded gene table from file: %s  (%d rows)\n', ...
                geneTable, height(geneTable));
    end

    % Validate column names
    if ~ismember(fromCol, geneTable.Properties.VariableNames)
        error('renameModelGenes:badColumn', ...
              'Column "%s" not found in geneTable.\nAvailable columns: %s', ...
              fromCol, strjoin(geneTable.Properties.VariableNames, ', '));
    end
    if ~ismember(toCol, geneTable.Properties.VariableNames)
        error('renameModelGenes:badColumn', ...
              'Column "%s" not found in geneTable.\nAvailable columns: %s', ...
              toCol, strjoin(geneTable.Properties.VariableNames, ', '));
    end

    % Build mapping: oldName → newName
    % Use only rows where both fromCol and toCol are non-empty. If a fromCol value maps to
    % multiple toCol values (isoforms), take the first non-empty one — gene names are
    % per-gene, not per-isoform.

    fromVals = geneTable.(fromCol);
    toVals   = geneTable.(toCol);

    % Convert to char cell if needed
    if ~iscell(fromVals), fromVals = cellstr(fromVals); end
    if ~iscell(toVals),   toVals   = cellstr(toVals);   end

    % Build map (first occurrence wins for duplicates). Only genes with a
    % non-empty toCol value are mapped; the rest are left unchanged.
    nameMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for i = 1:numel(fromVals)
        k = strtrim(fromVals{i});
        v = strtrim(toVals{i});
        if ~isempty(k) && ~isempty(v) && ~isKey(nameMap, k)
            nameMap(k) = v;
        end
    end

    % Standardize grRules before renaming to ensure consistent formatting for accurate replacement.
    fprintf('Removing unnecessary outer parentheses ...\n');
    model.grRules = removeUnnecessaryParentheses(model.grRules);
    
    fprintf('Standardizing grRules ...\n');
    [model.grRules, ~, indexes2check] = standardizeGrRules(model, true);
    if ~isempty(indexes2check)
        warning('renameModelGenes:potentialErrors', ...
            '%d grRules have potential issues. Review them before proceeding.', ...
            numel(indexes2check));
    end

    % Rename model.genes; genes without a mapping are left unchanged
    nGenes   = numel(model.genes);
    notFound = {};
    for i = 1:nGenes
        oldName = model.genes{i};
        if isKey(nameMap, oldName)
            model.genes{i} = nameMap(oldName);
        else
            notFound{end+1} = oldName; %#ok<AGROW>
        end
    end

    % Report genes that had no mapping and were left unchanged
    if ~isempty(notFound)
        warning('renameModelGenes:noMapping', ...
            '%d gene(s) in the model had no mapping in column "%s" and were left unchanged:\n  %s', ...
            numel(notFound), fromCol, strjoin(unique(notFound), ', '));
    end

    % Rename genes inside grRules Sort keys longest-first to avoid partial substitutions
    % (e.g. replace 'gene10' before 'gene1').
    allKeys    = keys(nameMap);
    allVals    = values(nameMap);
    [~, order] = sort(cellfun(@numel, allKeys), 'descend');
    allKeys    = allKeys(order);
    allVals    = allVals(order);
    
    nRxns = numel(model.grRules);
    for i = 1:nRxns
        rule = model.grRules{i};
        if isempty(rule)
            continue
        end
        for j = 1:numel(allKeys)
            oldG = allKeys{j};
            newG = allVals{j};
            % Use word-boundary replacement: only match oldG when surrounded
            % by non-word characters (spaces, parentheses, start/end of string)
            % to avoid replacing 'gene1' inside 'gene10'.
            rule = regexprep(rule, ['(?<![^\s()])' regexptranslate('escape', oldG) ...
                                    '(?![^\s()])'], newG);
        end
        model.grRules{i} = rule;
    end

    % Rebuild rxnGeneMat to match the updated gene list
    fprintf('Rebuilding rxnGeneMat ...\n');
    [model.grRules, model.rxnGeneMat, ~] = standardizeGrRules(model, true);
    
    fprintf('Done. %d genes renamed, %d left unchanged.\n', ...
            nGenes - numel(notFound), numel(notFound));
end

%--------------------------------------------------------------------------
%  Helper functions

function grRules = removeUnnecessaryParentheses(grRules)
% Remove unnecessary outer parentheses from grRules.
% e.g., '(A and (B or C))' → 'A and (B or C)'
%       '(A or B)' → 'A or B'
% but preserves '(A and B) or C' (outer parens are needed here)
%
% Checks if the opening parenthesis at position 1 closes before the end
% of the string. If yes, outer parens are necessary. If no, they can be
% removed recursively.

    nRxns = numel(grRules);
    for i = 1:nRxns
        rule = grRules{i};
        if isempty(rule)
            continue
        end
    
        % Recursively remove outer parentheses while they are unnecessary
        while length(rule) >= 2 && rule(1) == '(' && rule(end) == ')'
            % Count brackets to see when the opening paren closes
            parenCount = 0;
            closedEarly = false;
    
            for j = 1:length(rule)
                if rule(j) == '('
                    parenCount = parenCount + 1;
                elseif rule(j) == ')'
                    parenCount = parenCount - 1;
                    % If it reaches 0 before the end, outer parens are necessary
                    if parenCount == 0 && j < length(rule)
                        closedEarly = true;
                        break
                    end
                end
            end
    
            if closedEarly
                % Opening paren closes early → outer parens are needed
                break
            else
                % Opening paren only closes at the end → they're unnecessary
                rule = rule(2:end-1);
            end
        end
    
        grRules{i} = rule;
    end
end
