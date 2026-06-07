function [newModel, rxnToCheck]=expandModel(model)
% expandModel
%   Expands a model which uses several gene associations for one reaction.
%   Each such reaction is split into several reactions, each under the control
%   of only one gene.
%  
% Input:
%   model       model structure
% 
% Output:
%   newModel    model structure with separate reactions for iso-enzymes, where
%               the reaction ids are renamed as to id_EXP_1, id_EXP_2, etc. 
%   rxnToCheck  cell array with original reaction identifiers for those
%               that contained nested and/or-relationships in grRules.
%
%   NOTE: grRules strings that involve nested expressions of 'and' and 'or'
%         might not be parsed correctly if they are not standardized (if the
%         standardizeGrRules functions was not first run on the model). For
%         those reactions, it is therefore advisable to inspect the reactions in
%         rxnToCheck to confirm correct model expansion.
%
% Usage: [newModel, rxnToCheck]=expandModel(model)

%Check how many reactions we will create (the number of or:s in the GPRs).
%This way, we can preallocate all fields and save much computation time

numOrs = count(model.grRules, ' or ');
toAdd = sum(numOrs);
prevNumRxns = length(model.rxns);
rxnToCheck={};
if toAdd > 0
    %Calculate indices to expand
    %For example, if a reaction with index x has 2 or:s, meaning it has 3
    %reactions after the split, we should add two copies of this reaction
    %For fields that should just be copied to the new reactions, we just keep
    %track of that there are two copies, i.e., we add x x to this vector.
    %That is exactly what repelem does for us.
    cpyIndices = repelem(1:prevNumRxns, numOrs);
    
    %Copy all fields that should just be copied
    model.S=[model.S model.S(:,cpyIndices)];
    model.rxnNames=[model.rxnNames;model.rxnNames(cpyIndices)];
    model.lb=[model.lb;model.lb(cpyIndices)];
    model.ub=[model.ub;model.ub(cpyIndices)];
    model.rev=[model.rev;model.rev(cpyIndices)];
    model.c=[model.c;model.c(cpyIndices)];
    if isfield(model,'subSystems')
        model.subSystems=[model.subSystems;model.subSystems(cpyIndices)];
    end
    if isfield(model,'eccodes')
        model.eccodes=[model.eccodes;model.eccodes(cpyIndices)];
    end
    if isfield(model,'equations')
        model.equations=[model.equations;model.equations(cpyIndices)];
    end
    if isfield(model,'rxnMiriams')
        model.rxnMiriams=[model.rxnMiriams;model.rxnMiriams(cpyIndices)];
    end
    if isfield(model,'rxnComps')
        model.rxnComps=[model.rxnComps;model.rxnComps(cpyIndices)];
    end
    if isfield(model,'rxnFrom')
        model.rxnFrom=[model.rxnFrom;model.rxnFrom(cpyIndices)];
    end
    if isfield(model,'rxnNotes')
        model.rxnNotes=[model.rxnNotes;model.rxnNotes(cpyIndices)];
    end
    if isfield(model,'rxnReferences')
        model.rxnReferences=[model.rxnReferences;model.rxnReferences(cpyIndices)];
    end
    if isfield(model,'rxnConfidenceScores')
        model.rxnConfidenceScores=[model.rxnConfidenceScores;model.rxnConfidenceScores(cpyIndices)];
    end
    if isfield(model,'rxnDeltaG')
        model.rxnDeltaG=[model.rxnDeltaG;model.rxnDeltaG(cpyIndices)];
    end
    
    %now expand the more complex fields - will be filled in later
    model.rxns=[model.rxns;cell(toAdd,1)];
    model.grRules=[model.grRules;cell(toAdd,1)];
    model.rxnGeneMat=[model.rxnGeneMat;sparse(toAdd,size(model.rxnGeneMat,2))];
    
    %Loop throught those reactions and fill in the expanded data
    nextIndex = prevNumRxns + 1;
    for i=1:prevNumRxns
        if (numOrs(i) > 0)
            %Check that it doesn't contain nested 'and' and 'or' relations and
            %print a warning if it does
            if ~isempty(strfind(model.grRules{i},' and '))
                rxnToCheck{end+1,1}=model.rxns{i};
            end

            %Get rid of all '(' and ')' since I'm not looking at complex stuff
            %anyways
            geneString=model.grRules{i};
            geneString=strrep(geneString,'(','');
            geneString=strrep(geneString,')','');
            geneString=strrep(geneString,' or ',';');

            %Split the string into gene names
            geneNames=regexp(geneString,';','split');

            %Update the reaction to only use the first gene
            model.grRules{i}=['(' geneNames{1} ')'];
            %Find the gene in the gene list If ' and ' relationship, first
            %split the genes
            model.rxnGeneMat(i,:)=0;
            if ~isempty(strfind(geneNames(1),' and '))
                andGenes=regexp(geneNames{1},' and ','split');
                model.rxnGeneMat(i,ismember(model.genes,andGenes)) = 1;
            else
                [~, index]=ismember(geneNames(1),model.genes);
                model.rxnGeneMat(i,index)=1;
            end

            %Insert the reactions at the end of the model and without
            %allocating space. This is not nice, but ok for now
            for j=2:numel(geneNames)
                ind = nextIndex+j-2;
                model.rxns{ind}=[model.rxns{i} '_EXP_' num2str(j)];
                
                model.grRules{ind}=['(' geneNames{j} ')'];

                if ~isempty(strfind(geneNames(j),' and '))
                    andGenes=regexp(geneNames{j},' and ','split');
                    model.rxnGeneMat(ind,ismember(model.genes,andGenes)) = 1;
                else
                    model.rxnGeneMat(ind,ismember(model.genes,geneNames(j))) = 1;
                end
            end
            model.rxns{i}=[model.rxns{i}, '_EXP_1'];
            nextIndex = nextIndex + numOrs(i);
        end        
    end
    newModel=model;
else
    %There are no reactions to expand, return the model as is
    newModel=model;
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(newModel,true);
newModel.grRules = grRules;
newModel.rxnGeneMat = rxnGeneMat;
end
