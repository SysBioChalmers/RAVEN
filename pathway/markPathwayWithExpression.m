function pathway=markPathwayWithExpression(pathway,model,experiment,experimentOrder)
% markPathwayWithExpression
%   Adds gene expression data to a pathway object
%
%   pathway           map structure as generated by constructPathwayFromCelldesigner
%   model             a model structure
%   experiment        experiment structure as generated by getExpressionStructure
%   experimentOrder   cell array with two experiment names as defined in
%                     'experiment'
%
%   pathway           updated pathway object
%
% Usage: pathway=markPathwayWithExpression(pathway,model,experiment,experimentOrder)

if numel(experimentOrder)~=2
    EM='This can only be done for two cases at the moment (experimentOrder must be two elements)';
    dispEM(EM);
end

%Check that experiment fit with experimentOrder
[present, expIds]=ismember(experimentOrder,experiment.experiments);

if ~all(present)
    EM='Not all experiments could be found in the experiment structure';
    dispEM(EM);
end

experiment.data=experiment.data(:,expIds);
experiment.experiments=experiment.experiments(:,expIds);

%Go through each species in pathway and if it's a reaction, get it's genes
%and log2-fold change in expression from model and experiment
for i=1:numel(pathway.listOfSpecies)
    if strcmp('PROTEIN', pathway.listOfSpecies(i).type)
        if isfield(pathway.listOfSpecies(i),'note')
            if ~isempty(pathway.listOfSpecies(i).note)
                %Get the reaction if present in model
                [present, index]=ismember(pathway.listOfSpecies(i).note,model.rxns);
                
                %If present, then get the genes
                if any(present)
                    [~, genes]=find(model.rxnGeneMat(index,:));
                    
                    %If it was associated with genes match them to the ORFs
                    if any(genes)
                        [present, experimentIndexes]=ismember(model.genes(genes),experiment.orfs);
                        
                        %Add annotation to pathway structure
                        if any(present)
                            pathway.listOfSpecies(i).orfs=experiment.orfs(experimentIndexes(present));
                            pathway.listOfSpecies(i).expA=experiment.data(experimentIndexes(present),1);
                            pathway.listOfSpecies(i).expB=experiment.data(experimentIndexes(present),2);
                        end
                    end
                end
            end
        end
    end
end
end
