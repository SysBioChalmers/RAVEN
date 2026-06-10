function setTitle(handle, pathway, title)
% setTitle  Print a title on a map drawn by drawPathway.
%
% Parameters
% ----------
% handle : double
%     handle to the figure on which to print the title.
% pathway : struct
%     pathway structure representing the pathway to be drawn.
% title : char
%     the title.
%
% Examples
% --------
%     setTitle(handle, pathway, title);

dimension=getPathwayDimensions(pathway);

handle=text(dimension(1)+0.5*dimension(3), -150,...
    title,'fontname','Small Fonts','fontsize',5,...
    'interpreter', 'tex','HorizontalAlignment','center','verticalalignment','middle');
end
