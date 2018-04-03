function newModel=expandModel(model)
% expandModel
%   Expands a model which uses several gene associations for one reaction.
%   Each such reaction is split into several reactions, each under the control
%   of only one gene.
%
%   model     A model structure
%
%   newModel  A model structure with separate reactions for iso-enzymes
%
%	The reaction ids are renamed according to id_EXP_1, id_EXP_2..
%
%   NOTE: As it is now this code might not work for advanced grRules strings
%   that involve nested expressions of 'and' and 'or'
%
%   Usage: newModel=expandModel(model)
%
%   Simonas Marcisauskas, 2018-04-03
%

%Start by checking which reactions could be expanded
rxnsToExpand=false(numel(model.rxns),1);

for i=1:numel(model.rxns)
    if ~isempty(strfind(model.grRules{i},' or '))
        rxnsToExpand(i)=true;
    end
end

rxnsToExpand=find(rxnsToExpand);

if any(rxnsToExpand)
    %Loop throught those reactions and expand them
    for i=1:numel(rxnsToExpand)
        %Check that it doesn't contain nested 'and' and 'or' relations and
        %print a warning if it does
        if ~isempty(strfind(model.grRules{rxnsToExpand(i)},' and '))
            EM=['Reaction ' model.rxns{rxnsToExpand(i)} ' contains nested and/or-relations. Large risk of errors'];
            dispEM(EM,false);
        end

        %Get rid of all '(' and ')' since I'm not looking at complex stuff
        %anyways
        geneString=model.grRules{rxnsToExpand(i)};
        geneString=strrep(geneString,'(','');
        geneString=strrep(geneString,')','');
        geneString=strrep(geneString,' or ',';');

        %Split the string into gene names
        geneNames=regexp(geneString,';','split');

        %Update the reaction to only use the first gene
        model.grRules{rxnsToExpand(i)}=['(' geneNames{1} ')'];
        %Find the gene in the gene list
        %If ' and ' relationship, first split the genes
        model.rxnGeneMat(rxnsToExpand(i),:)=0;
        if ~isempty(strfind(geneNames(1),' and '))
            andGenes=regexp(geneNames{1},' and ','split');
            for h=1:numel(andGenes)
                [~, index]=ismember(andGenes(h),model.genes);
                model.rxnGeneMat(rxnsToExpand(i),index)=1;
            end
        else
            [~, index]=ismember(geneNames(1),model.genes);
            model.rxnGeneMat(rxnsToExpand(i),index)=1;
        end

        %Insert the reactions at the end of the model and without
        %allocating space. This is not nice, but ok for now
        for j=2:numel(geneNames)
            model.rxns=[model.rxns;[model.rxns{rxnsToExpand(i)} '_EXP_' num2str(j)]];
            model.rxnNames=[model.rxnNames;model.rxnNames(rxnsToExpand(i))];
            model.lb=[model.lb;model.lb(rxnsToExpand(i))];
            model.ub=[model.ub;model.ub(rxnsToExpand(i))];
            model.rev=[model.rev;model.rev(rxnsToExpand(i))];
            model.c=[model.c;model.c(rxnsToExpand(i))];
            model.S=[model.S model.S(:,rxnsToExpand(i))];
            model.grRules=[model.grRules;['(' geneNames{j} ')']];

            pad=zeros(1,numel(model.genes));
            if ~isempty(strfind(geneNames(j),' and '))
                andGenes=regexp(geneNames{j},' and ','split');
                for h=1:numel(andGenes)
                    [~, index]=ismember(andGenes(h),model.genes);
                    pad(index)=1;
                end
            else
                [~, index]=ismember(geneNames(j),model.genes);
                pad(index)=1;
            end
            model.rxnGeneMat=[model.rxnGeneMat;pad];

            if isfield(model,'subSystems')
                model.subSystems=[model.subSystems;model.subSystems(rxnsToExpand(i))];
            end
            if isfield(model,'eccodes')
                model.eccodes=[model.eccodes;model.eccodes(rxnsToExpand(i))];
            end
            if isfield(model,'equations')
                model.equations=[model.equations;model.equations(rxnsToExpand(i))];
            end
            if isfield(model,'rxnMiriams')
                model.rxnMiriams=[model.rxnMiriams;model.rxnMiriams(rxnsToExpand(i))];
            end
            if isfield(model,'rxnComps')
                model.rxnComps=[model.rxnComps;model.rxnComps(rxnsToExpand(i))];
            end
            if isfield(model,'rxnFrom')
                model.rxnFrom=[model.rxnFrom;model.rxnFrom(rxnsToExpand(i))];
            end
            if isfield(model,'rxnNotes')
                model.rxnNotes=[model.rxnNotes;model.rxnNotes(rxnsToExpand(i))];
            end
            if isfield(model,'rxnReferences')
                model.rxnReferences=[model.rxnReferences;model.rxnReferences(rxnsToExpand(i))];
            end
            if isfield(model,'rxnConfidenceScores')
                model.rxnConfidenceScores=[model.rxnConfidenceScores;model.rxnConfidenceScores(rxnsToExpand(i))];
            end
        end
    end
    newModel=model;
else
    %There are no reactions to expand, return the model as is
    newModel=model;
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(newModel);
newModel.grRules = grRules;
newModel.rxnGeneMat = rxnGeneMat;
end
