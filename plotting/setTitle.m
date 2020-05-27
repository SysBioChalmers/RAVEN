function setTitle(handle, pathway, title)
% plotLabels
%   Prints a title on on a map drawn by drawPathway.
%
%   handle      handle to the figure on which to print the labels
%   pathway     pathway structure representing the pathway to be drawn
%   title       the title
%
%   Usage: setTitle(handle, pathway, title)
%
%   Rasmus Agren, 2010-12-16
%

dimension=getPathwayDimensions(pathway);

handle=text(dimension(1)+0.5*dimension(3), -150,...
    title,'fontname','Small Fonts','fontsize',5,...
    'interpreter', 'tex','HorizontalAlignment','center','verticalalignment','middle');
end
