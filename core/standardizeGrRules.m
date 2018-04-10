function [grRules,rxnGeneMat,indexes2check] = standardizeGrRules(model,embeded)
% standardizeGrRules
%   Standardizes gene-rxn rules in a model according to the following
%       - No overall containing brackets
%       - Just enzyme complexes are enclosed into brackets
%       - ' and ' & ' or ' strings are strictly set to lowercases
%
%   A rxnGeneMat matrix consistent with the standardized grRules is created.
%
%   model        a model structure
%   embeded      TRUE if this function is called inside of another 
%                RAVEN function, default FALSE
%
%   grRules      [nRxns x 1] cell array with the standardized grRules
%   rxnGeneMat   [nRxns x nGenes]Sparse matrix consitent with the
%                standardized grRules
%   
%   If this function is going to be used in a model reconstruction or
%   modification pipeline it is recommended to run this function just
%   at the beginning of the process.
%
%   Usage: [grRules,rxnGeneMat,indexes2check]=standardizeGrRules(model,embeded)
%
%   Ivan Domenzain, 2018-04-10
%

%Preallocate fields
[~,n]      = size(model.S);
[g,~]      = size(model.genes);
rxnGeneMat = sparse(n,g);
grRules    = cell(n,1);

if nargin<2
    embeded = false;
end

if isfield(model,'grRules')
    %Search for potential logical errors in the grRules field
    indexes2check = findPotentialErrors(model.grRules,embeded);
    
    for i=1:length(model.grRules)
        originalSTR = model.grRules{i};
        %standardize the non-conflicting grRules
        if ~ismember(originalSTR,indexes2check)
            grRules{i}  = originalSTR;
            newSTR      = [];
            %Non-empty grRules are splitted in all their different isoenzymes
            genesSets = getSimpleGeneSets(originalSTR);

            if ~isempty(genesSets)
                for j=1:length(genesSets)
                    simpleSet  = genesSets{j};
                    rxnGeneMat = modifyRxnGeneMat(simpleSet,model.genes,rxnGeneMat,i);
                    %Enclose simpleSet in brackets
                    if length(genesSets)>1
                        if ~isempty(strfind(simpleSet,' and '))
                            simpleSet = horzcat('(',simpleSet,')');
                        end
                    end
                    %Separate genesSets in the substring (in case of
                    %isoenzymes)
                    if j<length(genesSets)
                        newSTR = [newSTR, simpleSet, ' or '];
                        %Add the last simpleSet
                    else
                        newSTR = [newSTR, simpleSet];
                    end

                end
                grRules{i} = newSTR;
            end
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that gets a cell array with all the simple geneSets in a given
%grRule string
function genesSets = getSimpleGeneSets(originalSTR)
genesSets  = [];
%If gene rule is not empty split in all its different isoenzymes
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
%Function that gets a simple genes set (single gene or enzyme complex)
%associated with the i-th reaction and modifies the correspondent row in the 
%rxnGeneMat accordingly.
function rxnGeneMat = modifyRxnGeneMat(genesSet,modelGenes,rxnGeneMat,i)
%Get individual genes
STR   = strrep(genesSet,') and (',' and ');
genes = strsplit(STR,' ');
for k=1:length(genes)
    if ~strcmpi(genes(k),' and ')
        genePos = find(strcmpi(modelGenes,genes(k)));
        if ~isempty(genePos)
            rxnGeneMat(i,genePos) = 1;
            %else
            %In this case the gene should be added to the
            %genes field (and to all of its dependencies)
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that gets the model field grRules and returns the indexes of the
%rules in which the pattern ") and (" is present. 
function indexes2check = findPotentialErrors(grRules,embeded)
indexes_l     = find(~cellfun(@isempty,strfind(grRules,') and (')));
indexes_U     = find(~cellfun(@isempty,strfind(grRules,') AND (')));
indexes2check = union(indexes_l,indexes_U);

if ~isempty(indexes2check)
    STR = ['") AND (" relationships found in' '\n' '\n'];
    for i=1:length(indexes2check)
        index = indexes2check(i);
        STR = [STR '  - grRule #' num2str(index) ': ' grRules{index} '\n'];
    end
    STR = [STR '\n'];
    if embeded
        STR = [STR,' If these grRules have been already curated then please'];
        STR = [STR,' ignore this message, otherwise run standardizeGrRules'];
        STR = [STR,' in the next way for getting further instructions:\n'];
        STR = [STR,'\n [~,~,indexes2check] = standardizeGrRules(model)\n'];
        warning(sprintf(STR))
    else
        STR = [STR,' This kind of relationship should only be present in'];
        STR = [STR,' reactions catalysed by complexes of isoenzymes e.g.\n\n'];
        STR = [STR,'  - (G1 or G2) and (G3 or G4)\n\n If this is the case'];
        STR = [STR,' please modify those grRules manually, writing all the'];
        STR = [STR,'  possible combinations e.g.\n\n'];
        STR = [STR,'  - (G1 and G3) or (G1 and G4) or (G2 and G3) or (G2 and G4)\n\n'];
        STR = [STR,' For other cases please modify the correspondent grRules'];
        STR = [STR,' avoiding:\n'];
        STR = [STR,'  1) Overall container brackets\n'];
        STR = [STR,'  2) Single unit enzymes enclosed into brackets\n'];
        STR = [STR,'  3) The usage of uppercases out of the gene names','\n'];
        STR = [STR,'  4) Unbalanced brackets','\n\n'];
        warning(sprintf(STR))
    end
end
end