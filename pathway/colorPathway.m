function returnPathway = colorPathway(pathway, reactionIDs, fluxes, referenceFluxes, ...
    cutOff, defaultColor, upColor, downColor)
% colorPathway
%	Calculates the color for those reactions that should be marked in the
%	map and returns an updated pathway structure.
%
%   pathway         pathway structure of the metabolic network
%   reactionsIDs    cell array with the names of the reactions in the model
%   fluxes          vector with flux values from the simulation
%   referenceFluxes vector with fluxes for the reference simulation
%   cutOff          a reaction is only colored if the absolute value of 
%                   at least one of the fluxes is above the cutoff value
%                   (opt, default 0)
%   defaultColor    a color in Matlab format to be used if there are no
%                   changes between the fluxes. This color is also used to 
%                   calculate the transition between the colors for up and
%                   down regulated fluxes (opt, default [1 1 1])
%   upColor         a color in Matlab format to be used if the flux is 
%                   larger than the reference flux (opt, default [0 1 0])
%   downColor       a color in Matlab format to be used if the flux is 
%                   smaller than the reference flux (opt, default [1 0 0])
%
%   returnPathway   updated pathway structure which contains coloring data
%       color       array indicating the coloring of the enzyme
%       signChange  true if the reaction has changed direction
%
%	Usage: returnPathway = colorPathway(pathway, reactionIDs, fluxes, referenceFluxes, ...
%           cutOff, defaultColor, upColor, downColor)
%
%   Rasmus Agren, 2010-12-16
%

if nargin<8
    downColor=[1 0 0];
end
if nargin<7
    upColor=[0 1 0];
end
if nargin<6
    defaultColor=[1 1 1];
end
if nargin<5
    cutOff=0;
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
                correspondingReaction=find(strcmpi(reactionIDs,pathway.listOfSpecies(i).note));
                if correspondingReaction>0
                    %Check if either flux is above the cutoff value
                    if abs(referenceFluxes(correspondingReaction))>=cutOff ||...
                            abs(fluxes(correspondingReaction))>=cutOff
                        %Calculate the corresponding color
                        [color, signChange]=getColorCodes(referenceFluxes(correspondingReaction)...
                            ,fluxes(correspondingReaction),1 , defaultColor, upColor, downColor);
                        if ~isempty(color)
                            returnPathway.listOfSpecies(i).color=color{1,1};
                            returnPathway.listOfSpecies(i).signChange=signChange{1,1};
                        end
                    end
                end
            end
        end
    end
end
end
