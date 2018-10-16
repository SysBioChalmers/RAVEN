function [grRules,rxnGeneMat,indexes2check] = standardizeGrRules(model,embedded)
% standardizeGrRules
%   Standardizes gene-rxn rules in a model according to the following
%       - No overall containing brackets
%       - Just enzyme complexes are enclosed into brackets
%       - ' and ' & ' or ' strings are strictly set to lowercases
%
%   A rxnGeneMat matrix consistent with the standardized grRules is created.
%
%   model        a model structure
%   embedded     true if this function is called inside of another 
%                RAVEN function (opt, default false)
%
%   grRules      [nRxns x 1] cell array with the standardized grRules
%   rxnGeneMat   [nRxns x nGenes]Sparse matrix consitent with the
%                standardized grRules
%   
%   If this function is going to be used in a model reconstruction or
%   modification pipeline it is recommended to run this function just
%   at the beginning of the process.
%
%   Usage: [grRules,rxnGeneMat,indexes2check]=standardizeGrRules(model,embedded)
%
%   Ivan Domenzain, 2018-05-04
%

%Preallocate fields
n          = length(model.rxns);
[g,~]      = size(model.genes);
rxnGeneMat = sparse(n,g);
grRules    = cell(n,1);
genes      = model.genes;
if nargin<2
    embedded = false;
end

if isfield(model,'grRules')
    originalGrRules=model.grRules; 
    originalGrRules=grRulesPreparation(originalGrRules);
    %Search for potential logical errors in the grRules field
    indexes2check = findPotentialErrors(originalGrRules,embedded,model);
    
    for i=1:length(originalGrRules)
        originalSTR = originalGrRules{i};
        grRules{i}  = originalSTR;
        %Non-empty grRules are splitted in all their different isoenzymes
        genesSets   = getSimpleGeneSets(originalSTR);
        rxnGeneMat  = modifyRxnGeneMat(genesSets,genes,rxnGeneMat,i);
        %standardize the non-conflicting grRules
        if ~ismember(i,indexes2check)
            newSTR = [];
            if ~isempty(genesSets)
                %For each simple genes set in the rule
                for j=1:length(genesSets)
                    simpleSet  = genesSets{j};
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
                %Update grRule
                grRules{i} = char(newSTR);
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
    %Remove all brackets
    originalSTR = strrep(originalSTR,'(','');
    originalSTR = strrep(originalSTR,')','');
    %Split all the different genesSets
    genesSets  = transpose(strsplit(originalSTR,' or '));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that gets a cell array of simple genes sets (single genes or
%enzyme complexes) associated with the i-th reaction and modifies the
%correspondent row in the rxnGeneMat accordingly.
function rxnGeneMat = modifyRxnGeneMat(genesSets,modelGenes,rxnGeneMat,i)

if ~isempty(genesSets)
    for j=1:length(genesSets)
        simpleSet  = genesSets{j};
        %        rxnGeneMat = modifyRxnGeneMat(simpleSet,model.genes,rxnGeneMat,i);
        %Get individual genes
        STR   = strrep(simpleSet,') and (',' and ');
        genes = strsplit(STR,' ');
        for k=1:length(genes)
            if ~strcmpi(genes(k),' and ')
                %Get gene index
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
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that gets the model field grRules and returns the indexes of the
%rules in which the pattern ") and (" is present.
function indexes2check = findPotentialErrors(grRules,embedded,model)
indxs_l       = find(~cellfun(@isempty,strfind(grRules,') and (')));
indxs_l_L     = find(~cellfun(@isempty,strfind(grRules,') and')));
indxs_l_R     = find(~cellfun(@isempty,strfind(grRules,'and (')));
indexes2check = vertcat(indxs_l,indxs_l_L,indxs_l_R);
indexes2check = unique(indexes2check);

if ~isempty(indexes2check)
    
    if embedded
        EM = 'Potentially problematic ") AND (" in the grRules for reaction(s): ';
        dispEM(EM,false,model.rxns(indexes2check),true)
    else
        STR = 'Potentially problematic ") AND (", ") AND" or "AND ("relat';
        STR = [STR,'ionships found in\n\n'];
        for i=1:length(indexes2check)
            index = indexes2check(i);
            STR = [STR '  - grRule #' num2str(index) ': ' grRules{index} '\n'];
        end
        STR = [STR,'\n This kind of relationships should only be present '];
        STR = [STR,'in  reactions catalysed by complexes of isoenzymes e'];
        STR = [STR,'.g.\n\n  - (G1 or G2) and (G3 or G4)\n\n For these c'];
        STR = [STR,'ases modify the grRules manually, writing all the po'];
        STR = [STR,'ssible combinations e.g.\n\n  - (G1 and G3) or (G1 a'];
        STR = [STR,'nd G4) or (G2 and G3) or (G2 and G4)\n\n For other c'];
        STR = [STR,'ases modify the correspondent grRules avoiding:\n\n '];
        STR = [STR,' 1) Overall container brackets, e.g.\n        "(G1 a'];
        STR = [STR,'nd G2)" should be "G1 and G2"\n\n  2) Single unit en'];
        STR = [STR,'zymes enclosed into brackets, e.g.\n        "(G1)" s'];
        STR = [STR,'hould be "G1"\n\n  3) The use of uppercases for logi'];
        STR = [STR,'cal operators, e.g.\n        "G1 OR G2" should be "G'];
        STR = [STR,'1 or G2"\n\n  4) Unbalanced brackets, e.g.\n        '];
        STR = [STR,'"((G1 and G2) or G3" should be "(G1 and G2) or G3"\n'];
        warning(sprintf(STR))
    end
end
end

function grRules = grRulesPreparation(grRules)
%Remove unnecessary blanks
grRules=strrep(grRules,'  ',' ');
grRules=strrep(grRules,'( ','(');
grRules=strrep(grRules,' )',')');
% Make sure that AND and OR strings are in lowercase
grRules=strrep(grRules,' AND ',' and ');
grRules=strrep(grRules,' OR ',' or ');
end