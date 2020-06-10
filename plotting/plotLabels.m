function plotLabels(handle,pathway)
% plotLabels
%   Plots labels on the pathway map.
%
%   handle      handle to the figure on which to print the labels
%   pathway     pathway structure of the metabolic network
%
%   NOTE: Since there is not dedicated shape for labels the PHENOTYPE
%           shape is used. The name of the PHENOTYPE object is used.
%
%   Usage: plotLabels(handle,pathway)

for i=1:numel(pathway.listOfSpecies)
    if strcmp(pathway.listOfSpecies(i).type,'PHENOTYPE')
        handle=text(pathway.listOfSpecies(i).x+pathway.listOfSpecies(i).w/2, pathway.listOfSpecies(i).y,...
            pathway.listOfSpecies(i).name,'fontname','Small Fonts','fontsize',4,...
            'interpreter', 'tex','HorizontalAlignment','center','verticalalignment','middle');
    end
end
end
