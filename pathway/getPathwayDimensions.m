function dimensions=getPathwayDimensions(pathway)
% getPathwayDimensions 
%   Retrieves the dimensions of metabolic network in a pathway structure.
%   Returns the position of the upper left corner, width and height.
%
%   pathway         pathway structure representing the pathway to be drawn
%
%   dimension       a 1x4 vector with x and y for the upper left corner,
%                   height and width
%
%   Usage: dimensions=getPathwayDimensions(pathway)
%
%   Rasmus Agren, 2010-12-16
%

right=0;
left=inf;
top=inf;
bottom=0;

%Loops through the compartments to find the right and bottom border and the
%position of the upper left corner
for i=1:length(pathway.listOfCompartments)
    if pathway.listOfCompartments(1,i).x<left
        left=pathway.listOfCompartments(1,i).x;
    end
    if pathway.listOfCompartments(1,i).y<top
        top=pathway.listOfCompartments(1,i).y;
    end
    if (pathway.listOfCompartments(1,i).x+pathway.listOfCompartments(1,i).w)>right
        right=pathway.listOfCompartments(1,i).x+pathway.listOfCompartments(1,i).w;
    end
    if (pathway.listOfCompartments(1,i).y+pathway.listOfCompartments(1,i).h)>bottom
        bottom=pathway.listOfCompartments(1,i).y+pathway.listOfCompartments(1,i).h;
    end
end

%Loops through the species to find the object furthest to the right, left, top
%and bottom
for i=1:length(pathway.listOfSpecies)
    if pathway.listOfSpecies(1,i).x<left
        left=pathway.listOfSpecies(1,i).x;
    end
    if pathway.listOfSpecies(1,i).y<top
        top=pathway.listOfSpecies(1,i).y;
    end
    if (pathway.listOfSpecies(1,i).x+pathway.listOfSpecies(1,i).w)>right
        right=pathway.listOfSpecies(1,i).x+pathway.listOfSpecies(1,i).w;
    end
    if (pathway.listOfSpecies(1,i).y+pathway.listOfSpecies(1,i).h)>bottom
        bottom=pathway.listOfSpecies(1,i).y+pathway.listOfSpecies(1,i).h;
    end
end

dimensions=[left,top,right-left,bottom-top];
end