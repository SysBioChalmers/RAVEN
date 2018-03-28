function newModel = standardizeGeneRules(model)
% standardizeGeneRules
%   Standardizes gene-rxn rules in a model and modifies rxnGeneMat for 
%   providing consistency with the grRules field
%
%   model        a model structure
%
%   newModel     an updated model structure
%
%   Usage: newModel = standardizeGeneRules(model)
%
%   Ivan Domenzain, 2018-03-28
%
newModel = model;
% Preallocate rxnGeneMat
[~,n]    = size(model.S);
[g,~]    = size(model.genes);
RGMat    = sparse(n,g);

if isfield(model,'grRules')
    % Search logical errors in the grRules field
    findLogicalErrors(model)  
    for i=1:length(model.grRules)
        originalSTR = model.grRules{i};
        newSTR      = [];
        % Non-empty grRules are splitted in all their different isoenzymes
        genesSets = getSimpleGeneSets(originalSTR);
        if ~isempty(genesSets)
            for j=1:length(genesSets)
                simpleSet  = genesSets{j};
                RGMat = modifyRxnGeneMat(simpleSet,model.genes,RGMat,i);
                % Enclose simpleSet in brackets
                if length(genesSets)>1
                    if ~isempty(strfind(simpleSet,' and '))
                        simpleSet = horzcat('(',simpleSet,')');
                    end
                end
                % Separate genesSets in the substring (in case of
                % isoenzymes)
                if j<length(genesSets)
                    newSTR = [newSTR, simpleSet, ' or '];
                    % Add the last simpleSet
                else
                    newSTR = [newSTR, simpleSet];
                end
                
            end
            newModel.grRules{i} = newSTR;
        end
    end
end
newModel.rxnGeneMat = RGMat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets a cell array with all the simple geneSets in a given 
% grRule string 
function genesSets = getSimpleGeneSets(originalSTR)
    genesSets  = [];
    % If gene rule is not empty split in all its different isoenzymes
    if ~isempty(originalSTR)
        originalSTR = strtrim(originalSTR);
        originalSTR = strrep(originalSTR,' OR ',' or ');
        originalSTR = strrep(originalSTR,' AND ',' and ');
        %Remove all brackets
        originalSTR = strrep(originalSTR,'(','');
        originalSTR = strrep(originalSTR,')','');
        %Split all the different genesSets
        genesSets  = transpose(strsplit(originalSTR,' or '));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets a simple genes set (single gene or enzyme complex) and
% a rxn index. The rxnGeneMat row (rxn) will be modified accordingly to the
% provided genes
function rxnGeneMat = modifyRxnGeneMat(genesSet,modelGenes,rxnGeneMat,i)
    %Get individual genes
    genes = strsplit(genesSet,' ');
    for k=1:length(genes)
        if ~strcmpi(genes(k),' and ')
            genePos = find(strcmpi(modelGenes,genes(k)));
            if ~isempty(genePos)
                rxnGeneMat(i,genePos) = 1;
                %else
                % In this case the gene should be added to the
                % genes field (and to all of its dependencies)
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets a simple genes set (single gene or enzyme complex) and
% a rxn index. The rxnGeneMat row (rxn) will be modified accordingly to the
function findLogicalErrors(model)  
    errors_l = find(~cellfun(@isempty,strfind(model.grRules,') and (')));
    errors_U = find(~cellfun(@isempty,strfind(model.grRules,') AND (')));
    errors   = union(errors_l,errors_U);
    if ~isempty(errors)
        for i=1:length(errors)
            index = errors(i);
            disp(['     grRule #:',num2str(index),' ',model.grRules{index}])
        end
        error('Logical errors found on grRules') 
    end
end