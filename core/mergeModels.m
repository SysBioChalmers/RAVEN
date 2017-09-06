function model=mergeModels(models,supressWarnings)
% mergeModels
%   Merges models into one model structure.
%
%   models          a cell array with model structures
%   supressWarnings true if warnings should be supressed (opt, default
%                   false)
%
%   model     a model structure with the merged model. Follows the structure
%             of normal models but also has 'rxnFrom/metFrom/geneFrom' fields
%             to indicate from which model each reaction/metabolite/gene was
%             taken
%
%   Usage: model=mergeModels(models)
%
%   Simonas Marcisauskas, 2017-08-25
%

%Just return the model
if numel(models)<=1
    model=models{1};
    return;
end

if nargin<2
    supressWarnings=false;
end

%Add new functionality in the order specified in models
model=models{1};
model.id='MERGED';
model.description='';

model.rxnFrom=cell(numel(models{1}.rxns),1);
model.rxnFrom(:)={models{1}.id};
model.metFrom=cell(numel(models{1}.mets),1);
model.metFrom(:)={models{1}.id};
if isfield(models{1},'genes')
    model.geneFrom=cell(numel(models{1}.genes),1);
    model.geneFrom(:)={models{1}.id};
end

if isfield(model,'subSystems')
    hasDeletedSubSystem=false;
else
    if supressWarnings==false
        EM='Cannot add subsystems since the existing model has no subsystems info. All reactions must have a subsystem for this to be included';
        dispEM(EM,false);
    end
    hasDeletedSubSystem=true;
end

if isfield(model,'equations')
    model=rmfield(model,'equations');
end

for i=2:numel(models)
    %Add the model id to the rxn id id it already exists in the model (id
    %have to be unique)
    %This is because it makes a '[]' string if no new reactions
    if ~isempty(models{i}.rxns)
        I=ismember(models{i}.rxns,model.rxns);
        models{i}.rxns(I)=strcat(models{i}.rxns(I),['_' models{i}.id]);
    end

    %Make sure that there are no conflicting reaction ids
    [~, ~, conflicting]=intersect(model.rxns,models{i}.rxns);

    if ~isempty(conflicting)
       printString=cell(numel(conflicting),1);
       for j=1:numel(conflicting)
           printString{j}=['Old: ' models{i}.rxns{conflicting(j)} ' New: ' models{i}.rxns{conflicting(j)} '_' models{i}.id];
           models{i}.rxns{conflicting(j)}=[models{i}.rxns{conflicting(j)} '_' models{i}.id];
       end
       if supressWarnings==false
           EM=['The following reaction IDs in ' models{i}.id ' are already present in the model and were renamed:'];
           dispEM(EM,false,printString);
           fprintf('\n');
       end
    end

    %Add all static stuff
    rxnFrom=cell(numel(models{i}.rxns),1);
    rxnFrom(:)={models{i}.id};
    model.rxnFrom=[model.rxnFrom;rxnFrom];
    model.rxns=[model.rxns;models{i}.rxns];
    model.rxnNames=[model.rxnNames;models{i}.rxnNames];
    model.lb=[model.lb;models{i}.lb];
    model.ub=[model.ub;models{i}.ub];
    model.c=[model.c;models{i}.c];
    model.rev=[model.rev;models{i}.rev];

    if hasDeletedSubSystem==false
        if isfield(models{i},'subSystems')
        	model.subSystems=[model.subSystems;models{i}.subSystems];
        else
           if supressWarnings==false
               EM='Cannot add subsystems since the existing model has no subsystems info. All reactions must have a subsystem for this to be included. Deleting subSystems field';
               dispEM(EM,false);
           end
           hasDeletedSubSystem=true;
           model=rmfield(model,'subSystems');
        end
    end

    if isfield(models{i},'eccodes')
       if isfield(model,'eccodes')
           model.eccodes=[model.eccodes;models{i}.eccodes];
       else
           emptyEC=cell(numel(model.rxns)-numel(models{i}.rxns),1);
           emptyEC(:)={''};
           model.eccodes=[emptyEC;models{i}.eccodes];
       end
    else
       if isfield(model,'eccodes')
           emptyEC=cell(numel(models{i}.rxns),1);
           emptyEC(:)={''};
           model.eccodes=[model.eccodes;emptyEC];
       end
    end

    if isfield(models{i},'rxnMiriams')
       if isfield(model,'rxnMiriams')
           model.rxnMiriams=[model.rxnMiriams;models{i}.rxnMiriams];
       else
           model.rxnMiriams=[cell(numel(model.rxns)-numel(models{i}.rxns),1);models{i}.rxnMiriams];
       end
    else
       if isfield(model,'rxnMiriams')
           model.rxnMiriams=[model.rxnMiriams;cell(numel(models{i}.rxns),1)];
       end
    end

    if isfield(models{i},'rxnNotes')
       if isfield(model,'rxnNotes')
           model.rxnNotes=[model.rxnNotes;models{i}.rxnNotes];
       else
           emptyNotes=cell(numel(model.rxns)-numel(models{i}.rxns),1);
           emptyNotes(:)={''};
           model.rxnNotes=[emptyNotes;models{i}.rxnNotes];
       end
    else
       if isfield(model,'rxnNotes')
           emptyNotes=cell(numel(models{i}.rxns),1);
           emptyNotes(:)={''};
           model.rxnNotes=[model.rxnNotes;emptyNotes];
       end
    end

    if isfield(models{i},'rxnReferences')
       if isfield(model,'rxnReferences')
           model.rxnReferences=[model.rxnReferences;models{i}.rxnReferences];
       else
           emptyReferences=cell(numel(model.rxns)-numel(models{i}.rxns),1);
           emptyReferences(:)={''};
           model.rxnReferences=[emptyReferences;models{i}.rxnReferences];
       end
    else
       if isfield(model,'rxnReferences')
           emptyReferences=cell(numel(models{i}.rxns),1);
           emptyReferences(:)={''};
           model.rxnReferences=[model.rxnReferences;emptyReferences];
       end
    end

    if isfield(models{i},'rxnConfidenceScores')
       if isfield(model,'rxnConfidenceScores')
           model.rxnConfidenceScores=[model.rxnConfidenceScores;models{i}.rxnConfidenceScores];
       else
           emptyConfidenceScores=cell(numel(model.rxns)-numel(models{i}.rxns),1);
           emptyConfidenceScores(:)={''};
           model.rxnConfidenceScores=[emptyConfidenceScores;models{i}.rxnConfidenceScores];
       end
    else
       if isfield(model,'rxnConfidenceScores')
           emptyConfidenceScores=cell(numel(models{i}.rxns),1);
           emptyConfidenceScores(:)={''};
           model.rxnConfidenceScores=[model.rxnConfidenceScores;emptyConfidenceScores];
       end
    end

    if isfield(models{i},'rxnComps')
       if isfield(model,'rxnComps')
           model.rxnComps=[model.rxnComps;models{i}.rxnComps];
       else
           model.rxnComps=[ones(numel(model.rxns)-numel(models{i}.rxns),1);models{i}.rxnComps];
           fprintf('NOTE: One of the models does not contain compartment information for its reactions. All reactions in that model has been assigned to the first compartment\n');
       end
    else
       if isfield(model,'rxnComps')
           model.rxnComps=[model.rxnComps;ones(numel(models{i}.rxns),1)];
           fprintf('NOTE: One of the models does not contain compartment information for its reactions. All reactions in that model has been assigned to the first compartment\n');
       end
    end

    if isfield(models{i},'rxnScores')
       if isfield(model,'rxnScores')
           model.rxnScores=[model.rxnScores;models{i}.rxnScores];
       else
           emptyRS=zeros(numel(model.rxns)-numel(models{i}.rxns),1);
           model.rxnScores=[emptyRS;models{i}.rxnScores];
       end
    else
       if isfield(model,'rxnScores')
           emptyRS=zeros(numel(models{i}.rxns),1);
           model.rxnScores=[model.rxnScores;emptyRS];
       end
    end

    %Get the new metabolites from matching the models.
    %Metabolites are said to be the same if they share name and
    %compartment id. This means that metabolite IDs are not taken into
    %account.
    oldMetComps=model.comps(model.metComps);
    oldMets=strcat(model.metNames,'[',oldMetComps,']');
    %This is because it makes a '[]' string if no new metabolites
    if ~isempty(models{i}.metNames)
        newMetComps=models{i}.comps(models{i}.metComps);
        newMets=strcat(models{i}.metNames,'[',newMetComps,']');
    else
        newMets={};
        newMetComps={};
    end
    tf=ismember(newMets,oldMets);
    metsToAdd=find(~tf);

    %First add the new metabolites
    %Make sure that there are no conflicting metabolite ids
    conflicting=ismember(models{i}.mets(metsToAdd),model.mets);

    conflicting=find(conflicting);

    if ~isempty(conflicting)
       printString=cell(numel(conflicting),1);
       for j=1:numel(conflicting)
           printString{j}=['Old: ' models{i}.mets{metsToAdd(conflicting(j))} ' New: ' models{i}.mets{metsToAdd(conflicting(j))} '_' models{i}.id];
           models{i}.mets{metsToAdd(conflicting(j))}=[models{i}.mets{metsToAdd(conflicting(j))} '_' models{i}.id];
       end
       if supressWarnings==false
           EM=['The following metabolite IDs in ' models{i}.id ' are already present in the model and were renamed:'];
           dispEM(EM,false,printString);
       end
    end

    %Add static info on the metabolites
    metFrom=cell(numel(metsToAdd),1);
    metFrom(:)={models{i}.id};
    model.metFrom=[model.metFrom;metFrom];
    model.mets=[model.mets;models{i}.mets(metsToAdd)];
    model.metNames=[model.metNames;models{i}.metNames(metsToAdd)];
    model.b=[model.b;zeros(numel(metsToAdd),size(model.b,2))];

    if isfield(model,'unconstrained')
       if isfield(models{i},'unconstrained')
            model.unconstrained=[model.unconstrained;models{i}.unconstrained(metsToAdd)];
       else
           model.unconstrained=[model.unconstrained;zeros(numel(metsToAdd),1)];
       end
    else
       if isfield(models{i},'unconstrained')
          model.unconstrained=[zeros(numel(model.mets),1);models{i}.unconstrained(metsToAdd)];
       end
    end

    %Only add extra info on new metabolites since it's a little tricky to
    %chose what to keep otherwise. Should change in the future
    if ~isempty(metsToAdd)
        if isfield(models{i},'inchis')
           if isfield(model,'inchis')
               model.inchis=[model.inchis;models{i}.inchis(metsToAdd)];
           else
               emptyInchi=cell(numel(model.mets)-numel(metsToAdd),1);
               emptyInchi(:)={''};
               model.inchis=[emptyInchi;models{i}.inchis(metsToAdd)];
           end
        else
           if isfield(model,'inchis')
               emptyInchi=cell(numel(metsToAdd),1);
               emptyInchi(:)={''};
               model.inchis=[model.inchis;emptyInchi];
           end
        end

        if isfield(models{i},'metFormulas')
           if isfield(model,'metFormulas')
               model.metFormulas=[model.metFormulas;models{i}.metFormulas(metsToAdd)];
           else
               emptyMetFormulas=cell(numel(model.mets)-numel(metsToAdd),1);
               emptyMetFormulas(:)={''};
               model.metFormulas=[emptyMetFormulas;models{i}.metFormulas(metsToAdd)];
           end
        else
           if isfield(model,'metFormulas')
               emptyMetFormulas=cell(numel(metsToAdd),1);
               emptyMetFormulas(:)={''};
               model.metFormulas=[model.metFormulas;emptyMetFormulas];
           end
        end

        if isfield(models{i},'metCharges')
           if isfield(model,'metCharges')
               model.metCharges=[model.metCharges;models{i}.metCharges(metsToAdd)];
           else
               emptyMetCharge=cell(numel(model.mets)-numel(metsToAdd),1);
               emptyMetCharge(:)={''};
               model.metCharges=[emptyMetCharge;models{i}.metCharges(metsToAdd)];
           end
        else
           if isfield(model,'metCharges')
               emptyMetCharge=cell(numel(metsToAdd),1);
               emptyMetCharge(:)={''};
               model.metCharges=[model.metCharges;emptyMetCharge];
           end
        end

        if isfield(models{i},'metMiriams')
           if isfield(model,'metMiriams')
               model.metMiriams=[model.metMiriams;models{i}.metMiriams(metsToAdd)];
           else
               emptyMetMiriam=cell(numel(model.mets)-numel(metsToAdd),1);
               model.metMiriams=[emptyMetMiriam;models{i}.metMiriams(metsToAdd)];
           end
        else
           if isfield(model,'metMiriams')
               emptyMetMiriam=cell(numel(metsToAdd),1);
               model.metMiriams=[model.metMiriams;emptyMetMiriam];
           end
        end
    end

    %Add if there are any new compartments and add those. This can change
    %the order of compartments and the corresponding indexes in
    %model.metComps.

    %Find overlapping and new compartments
    [overlap, oldIDs]=ismember(models{i}.comps,model.comps);
    overlap=find(overlap);

    %Add the new compartments if any
    if numel(overlap)~=numel(models{i}.compNames)
        compIndexes=oldIDs==0;

        %Make sure that there are no conflicting compartment ids
        [~, conflicting]=ismember(models{i}.compNames(compIndexes),model.compNames);
        EM=['The following compartment IDs in ' models{i}.id ' are already present in the model but with another name. They have to be renamed'];
        dispEM(EM,true,model.comps(conflicting));

        %It's ok to add duplicate name, but not duplicate IDs
        model.compNames=[model.compNames; models{i}.compNames(compIndexes)];
        model.comps=[model.comps; models{i}.comps(compIndexes)];
        if isfield(model,'compOutside')
            if isfield(models{i},'compOutside')
                model.compOutside=[model.compOutside; models{i}.compOutside(compIndexes)];
            else
                %This is if not all models have the field
                padding=cell(sum(compIndexes),1);
                padding(:)={''};
                model.compOutside=[model.compOutside;padding];
            end
        end
        if isfield(model,'compMiriams')
            if isfield(models{i},'compMiriams')
                model.compMiriams=[model.compMiriams; models{i}.compMiriams(compIndexes)];
            else
                %This is if not all models have the field
                model.compMiriams=[model.compMiriams;cell(sum(compIndexes),1)];
            end
        end
    end

    %Only add new comp info on the un-matched metabolites since the old ones will
    %be mapped to the existing list anyways
    [I, J]=ismember(newMetComps(metsToAdd),model.comps);
    %Just a check
    if ~all(I)
        EM='There was an unexpected error in matching compartments';
        dispEM(EM);
    end
    model.metComps=[model.metComps;J];

    %Create the new stoichiometric matrix
    model.S=[model.S;sparse(numel(metsToAdd),size(model.S,2))];

    %Rematch metabolite names. Not the most clever way to do it maybe
    allMets=strcat(model.metNames,'[',model.comps(model.metComps),']');
    [~, J]=ismember(newMets,allMets);

    %Update the stoichiometric matrix for the model to add
    newS=sparse(numel(model.mets),numel(models{i}.rxns));
    newS(J,:)=models{i}.S;

    model.S=[model.S newS];

    %Now add new genes
    if isfield(models{i},'genes')
        if ~isfield(model,'genes')
            %If there was no gene info before
            model.genes=models{i}.genes;
            model.rxnGeneMat=[sparse(numel(model.rxns),numel(models{i}.genes));models{i}.rxnGeneMat];
            emptyGene=cell(numel(model.rxns),1);
            emptyGene(:)={''};
            model.grRules=[emptyGene;models{i}.grRules];
            model.geneFrom=cell(numel(models{i}.genes),1);
            model.geneFrom(:)={models{i}.id};

            if isfield(models{i},'geneShortNames')
               model.geneShortNames=models{i}.geneShortNames;
            end

            if isfield(models{i},'geneMiriams')
               model.geneMiriams=models{i}.geneMiriams;
            end

            if isfield(models{i},'geneComps')
               model.geneComps=models{i}.geneComps;
            end
        else
            %If gene info should be merged
            a=ismember(models{i}.genes,model.genes);

            genesToAdd=find(~a);

            %Only add extra gene info on new genes. This might not be
            %correct and should be changed later...
            if ~isempty(genesToAdd)
                model.genes=[model.genes;models{i}.genes(genesToAdd)];
                emptyGene=cell(numel(genesToAdd),1);
                emptyGene(:)={models{i}.id};
                model.geneFrom=[model.geneFrom;emptyGene];
                model.rxnGeneMat=[model.rxnGeneMat sparse(size(model.rxnGeneMat,1),numel(genesToAdd))];

                if isfield(models{i},'geneShortNames')
                    if isfield(model,'geneShortNames')
                        model.geneShortNames=[model.geneShortNames;models{i}.geneShortNames(genesToAdd)];
                    else
                        emptyGeneSN=cell(numel(model.genes)-numel(genesToAdd),1);
                        emptyGeneSN(:)={''};
                        model.geneShortNames=[emptyGeneSN;models{i}.geneShortNames(genesToAdd)];
                    end
                else
                    if isfield(model,'geneShortNames')
                        emptyGeneSN=cell(numel(genesToAdd),1);
                        emptyGeneSN(:)={''};
                        model.geneShortNames=[model.geneShortNames;emptyGeneSN];
                    end
                end

                if isfield(models{i},'geneMiriams')
                    if isfield(model,'geneMiriams')
                        model.geneMiriams=[model.geneMiriams;models{i}.geneMiriams(genesToAdd)];
                    else
                        emptyGeneMir=cell(numel(model.genes)-numel(genesToAdd),1);
                        model.geneMiriams=[emptyGeneMir;models{i}.geneMiriams(genesToAdd)];
                    end
                else
                    if isfield(model,'geneMiriams')
                        emptyGeneMir=cell(numel(genesToAdd),1);
                        model.geneMiriams=[model.geneMiriams;emptyGeneMir];
                    end
                end

                if isfield(models{i},'geneComps')
                    if isfield(model,'geneComps')
                        model.geneComps=[model.geneComps;models{i}.geneComps(genesToAdd)];
                    else
                        emptyGeneMir=ones(numel(model.genes)-numel(genesToAdd),1);
                        model.geneComps=[emptyGeneMir;models{i}.geneComps(genesToAdd)];
                        EM='Adding genes with compartment information to a model without such information. All existing genes will be assigned to the first compartment';
                        dispEM(EM,false);
                    end
                else
                    if isfield(model,'geneComps')
                        emptyGeneMir=ones(numel(genesToAdd),1);
                        model.geneComps=[model.geneComps;emptyGeneMir];
                        EM='Adding genes with compartment information to a model without such information. All existing genes will be assigned to the first compartment';
                        dispEM(EM,false);
                    end
                end
            end

            %Remap the genes from the new model. The same thing as with
            %mets; this is a wasteful way to do it but I don't care right
            %now
            [a, b]=ismember(models{i}.genes,model.genes);

            %Just a check
            if ~all(a)
                EM='There was an unexpected error in matching genes';
                dispEM(EM);
            end

            %Create the new rxnGene matrix
            rxnGeneMat=sparse(numel(models{i}.rxns),numel(model.genes));
            rxnGeneMat(:,b)=models{i}.rxnGeneMat;
            model.rxnGeneMat=[model.rxnGeneMat; rxnGeneMat];
            model.grRules=[model.grRules;models{i}.grRules];
        end
    else
        %Add empty gene associations
        if isfield(model,'genes')
            model.rxnGeneMat=[model.rxnGeneMat;sparse(numel(models{i}.rxns),numel(model.genes))];
            emptyGene=cell(numel(models{i}.rxns),1);
            emptyGene(:)={''};
            model.grRules=[model.grRules;emptyGene];
        end
    end
end
end
