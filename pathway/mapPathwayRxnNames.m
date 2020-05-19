function [pathway, notMapped]=mapPathwayRxnNames(pathway,originalLabels,newLabels)
% mapPathwayRxnNames
%   For mapping labels in the pathway object. Useful if you want to change
%   what is shown in the reaction boxes.
%
%   pathway         pathway structure representing the pathway to be drawn
%   originalLabels  cell array with the original reaction labels
%   newLabels       cell array with the new reaction labels
%
%   pathway         an updated pathway structure
%   notMapped       a cell array with labels that could not be found in the
%                   pathway object
%
%   Usage: [pathway, notMapped]=mapPathwayRxnNames(pathway,originalLabels,newLabels)
%
%   Rasmus Agren, 2014-01-09
%

if numel(originalLabels)~=numel(newLabels)
    EM='The new label cell array must have the same length as the old label cell array';
    dispEM(EM);
end

mapped=false(numel(originalLabels),1);

for i=1:numel(pathway.listOfSpecies)
    if strcmpi(pathway.listOfSpecies(i).type,'PROTEIN')
        I=find(ismember(originalLabels,pathway.listOfSpecies(i).name));
        if any(I)
            if numel(I)==1
                pathway.listOfSpecies(i).name=newLabels{I};
                mapped(I)=true;
            else
                EM=['The label "' pathway.listOfSpecies(i).name '" was found in several positions in oldLabels'];
                dispEM(EM);
            end
        end
    end
end

notMapped=originalLabels(~mapped);
end
