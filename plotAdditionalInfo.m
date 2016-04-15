function plotAdditionalInfo(handle, pathway, additionalText, exampleBoxText,...
    maxChange, defaultColor, upColor, downColor)
% plotAdditionalInfo
%	Plots some additional information in a figure.
%
%   handle              handle to the figure to plot the information on
%   pathway             pathway structure of the metabolic network
%   additionalText      array list with additional text to print (e.g.
%                       fluxes or constraints)
%   exampleBoxText      array list of text to put in an example box. The
%                       text should explain what is printed in the enzyme
%                       boxes
%   maxChange           the logfold increase or decrease that corresponds
%                       to full negative or full positive coloring. Must
%                       be a positive value (opt, default 1)
%   defaultColor        a color in Matlab format to be used if there are no
%                       changes between the fluxes. This color is also used to 
%                       calculate the transition between the colors for up and
%                       down regulated fluxes (opt, default [1 1 1])
%   upColor             a color in Matlab format to be used if the flux is 
%                       larger than the reference flux (opt, default [0 1 0])
%   downColor           a color in Matlab format to be used if the flux is 
%                       smaller than the reference flux (opt, default [1 0
%                       0])
%
%   NOTE: At the moment the positions of the text/figures are (semi-)hard
%           coded.
%
%   Usage:  errorFlag = plotAdditionalInfo(handle, additionalText, exampleBoxText,...
%               maxChange, defaultColor, upColor, downColor)
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
    maxChange=1;
end

%Finds the the dimensions of the metabolic network. The additional information will be 
%positioned relative to that object.
dimension=getPathwayDimensions(pathway);

%Plots the example box
rectangle('edgecolor',[0 0 0], 'facecolor', defaultColor, 'linewidth', 1,...
    'position', [dimension(1)+dimension(3)+100 100 700 320],'curvature', [0.1 0.1]);
handle=text(dimension(1)+dimension(3)+100+6, 100+0.5*320,...
            exampleBoxText,'fontname','Small Fonts','fontsize',2,...
            'interpreter', 'tex','verticalalignment','middle','HorizontalAlignment','left');
handle=text(dimension(1)+dimension(3)+300+126, 20,...
            'EXAMPLE:','fontname','Small Fonts','fontsize',4,...
            'interpreter', 'tex','HorizontalAlignment','center','verticalalignment','middle');

%Calculates 10 colors between upColor and defaultColor
%The color is linear from the upColor to the defaultColor
colorValues=[];
for i=1:11
    logvalue=maxChange-(i-1)*0.1*maxChange;
	colorValues=[colorValues;...
                  [defaultColor(1)+(upColor(1)-defaultColor(1))*logvalue/(maxChange)...
                   defaultColor(2)+(upColor(2)-defaultColor(2))*logvalue/(maxChange)...
                   defaultColor(3)+(upColor(3)-defaultColor(3))*logvalue/(maxChange)]];
end
%The color is linear from the defaultColor to downColor
for i=1:10
    logvalue=i*0.1*maxChange;
	colorValues=[colorValues;...
                  [defaultColor(1)+(downColor(1)-defaultColor(1))*logvalue/(maxChange)...
                   defaultColor(2)+(downColor(2)-defaultColor(2))*logvalue/(maxChange)...
                   defaultColor(3)+(downColor(3)-defaultColor(3))*logvalue/(maxChange)]];
end

%Draw lines that represent the different colors
handle=text(dimension(1)+dimension(3)+100, 560,...
            'log10(|condition A|/|condition B|)','fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');
lengthLine=40;
startX=dimension(1)+dimension(3)+150;

startY=650;
width=6;
for i=1:size(colorValues,1)
    line([startX; startX],[startY+(i-1)*lengthLine; startY+i*lengthLine],'color',colorValues(i,:),'linewidth',width);
end

handle=text(startX+3*width, startY,...
            ['-  ' num2str(maxChange)],'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');
        
handle=text(startX+3*width, startY+5.5*lengthLine,...
            ['-  ' num2str(maxChange/2)],'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');

handle=text(startX+3*width, startY+10.5*lengthLine,...
            ['-  0'],'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');
        
handle=text(startX+3*width, startY+15.5*lengthLine,...
            ['-  -' num2str(maxChange/2)],'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');
        
handle=text(startX+3*width, startY+20.5*lengthLine,...
            ['-  -' num2str(maxChange)],'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'tex','HorizontalAlignment','left','verticalalignment','middle');

%Plots the additional text
handle=text(dimension(1)+dimension(3)+100, 1800,...
            additionalText,'fontname','Small Fonts','fontsize',3,...
            'interpreter', 'none','verticalalignment','top','HorizontalAlignment','left');