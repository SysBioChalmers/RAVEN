function [returnPathway, errorFlag]= markPathwayWithFluxes(pathway, reactionIDs, fluxes, referenceFluxes)
% markPathwayWithFluxes
%   Marks each enzyme in a pathway structure with the corresponding fluxes
%   from two simulation results. This is done for enzymes that has the name
%   of a reaction in the note field. The reaction has to be present in
%   reactionIDs.
%
%   pathway         pathway structure of the metabolic network
%   reactionsIDs    cell array with the names of the reactions in the model
%   fluxes          vector with flux values
%   referenceFluxes vector with fluxes to compare to
%
%   returnPathway   updates the original pathway structure by adding the
%                   fields flux and referenceFlux for each marked reaction
%   errorFlag       true if there has been an error
%
%   Usage: [returnPathway, errorFlag] = markPathwayWithFluxes(pathway, reactionIDs,
%   fluxes, referenceFluxes)

%Check if all variables are of the correct dimension
if length(reactionIDs)~=length(fluxes) || length(reactionIDs)~=length(referenceFluxes)
    returnPathway=pathway;
    errorFlag=1;
    fprintf('reactionIDs, fluxes, and referenceFluxes do not have the same number of elements.');
    return;
end

%Loop through the components in the pathway. Check if any component has a
%note and then check for the corresponding reaction id
returnPathway=pathway;

for i=1:length(pathway.listOfSpecies)
    if strcmpi(pathway.listOfSpecies(i).type,'PROTEIN')
        if isfield(pathway.listOfSpecies(i),'note')
            if ~isempty(pathway.listOfSpecies(i).note)
                %If there is a note check if there is a corresponding
                %reaction id
                index=find(strcmpi(reactionIDs,pathway.listOfSpecies(i).note));
                %If there is a match
                if any(index)
                    returnPathway.listOfSpecies(i).flux=fluxes(index);
                    returnPathway.listOfSpecies(i).referenceFlux=referenceFluxes(index);
                end
            end
        end
    end
end
end
