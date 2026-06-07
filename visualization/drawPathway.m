function drawPathway(pathway, h, cutOff, defaultColor)
% drawPathway
%	Draws a metabolic network to a figure.
%
%   pathway         pathway structure representing the pathway to be drawn
%   handle          handle to a figure (optional)
%   cutOff          the fluxes are only printed if the absolute value of
%                   at least one of the fluxes is above the cutoff value
%                   (optional, default 0)
%   defaultColor    color in Matlab format to be used as the background
%                   color for enzymes if no color is specified in the
%                   pathway structure (optional, default [1 1 1])
%
%   NOTE:   At the moment all text sizes and some positions are hard coded.
%           This means that this code is not appliable for any map
%
% Usage: drawPathway(pathway, handle, cutOff, defaultColor)

if nargin<4
    defaultColor=[1 1 1];
end
if nargin<3
    cutOff=0;
end
if nargin<2
    if ishandle(gcf)
        h=gcf;
    else
        h=figure();
    end
end

%Will contain all metabolite info and where to print them
metText={};
metX=[];
metY=[];
metRec=[];

%Same for reaction text (names, fluxes and maybe other)
rxnText={};
rxnX=[];
rxnY=[];
rxnLinesX=[];
rxnLinesY=[];
rxnRec=[];
rxnFace=[];
rxnEdge=[];
rxnArrows=[];
rxnArrowsColor=[];

%Expression boxes
expRec=[];
expFace=[];

%Lines to enzymes
enzLinesX=[];
enzLinesY=[];

%Same for compartments
compText={};
compX=[];
compY=[];
compRec=[];

%Check if it should plot expression data
if isfield(pathway.listOfSpecies(1),'orfs')
    drawExpression=true;
else
    drawExpression=false;
end

%Get the dimensions of the network
dimension=getPathwayDimensions(pathway);

%Specify the properties of the figure
set(h, 'PaperOrientation', 'landscape');
set(h, 'Color', [1 1 1]);
xmargin=0.83;
ymargin=1.32;
set(h, 'PaperPosition', [xmargin ymargin 29.67743169791-2*xmargin 20.98404194812-2*ymargin]);
set(h, 'Renderer', 'painters');
set(h, 'Position', [10 10 1000 700]);

%The height has to be at least the height of the information panel
%(approximately 3000)
height=dimension(2)+dimension(4)+200;
height=max(height,3000);

axes('ytick',[],'ydir','reverse','xtick',[],'ycolor',[1 1 1], 'xcolor',...
    [1 1 1],'position',[0 0 1 1], 'xLim',[dimension(1)-40 dimension(1)+dimension(3)+800],...
    'yLim', [dimension(2)-200 height]);
daspect([1,1,1]);

%Loops through the compartments.
for i=1:length(pathway.listOfCompartments)
    position=[pathway.listOfCompartments(1,i).x pathway.listOfCompartments(1,i).y...
        pathway.listOfCompartments(1,i).w pathway.listOfCompartments(1,i).h];
    compRec=[compRec;position];
    compText=[compText;upper(pathway.listOfCompartments(i).name)];
    compX=[compX;position(1)+0.5*position(3)];
    compY=[compY;position(2)+position(4)-70];
end

%Get positions of lines
for i=1:length(pathway.listOfReactions)
    for j=1:length(pathway.listOfReactions(i).componentList)
        %The line is always drawn from the middle to the component
        middle=pathway.listOfReactions(i).middlePoint;
        x(1)=middle(1);
        y(1)=middle(2);
        x(2)=pathway.listOfReactions(i).componentList(j).anchor(1);
        y(2)=pathway.listOfReactions(i).componentList(j).anchor(2);
        
        %Check to see if a line or an arrow should be drawn
        if strcmpi(pathway.listOfReactions(i).componentList(j).toArrow,'true')
            rxnArrows=[rxnArrows;x(1) y(1) x(2) y(2)];
            
            %Draw a red arrow if it's a base product
            if strcmpi(pathway.listOfReactions(i).componentList(j).baseProduct,'true')
                rxnArrowsColor=[rxnArrowsColor;1 0 0];
            else
                rxnArrowsColor=[rxnArrowsColor;0 0 0];
            end
        else
            %If the component is an enzyme then a different line should be
            %used
            if strcmpi(pathway.listOfReactions(i).componentList(j).type,'ENZYME')
                enzLinesX=[enzLinesX;x(1) x(2)];
                enzLinesY=[enzLinesY;y(1) y(2)];
            else
                rxnLinesX=[rxnLinesX;x(1) x(2)];
                rxnLinesY=[rxnLinesY;y(1) y(2)];
            end
        end
    end
end

%Loops through the species
for i=1:length(pathway.listOfSpecies)
    %Species should be represented with ellipses and enzymes with
    %rectangles
    position=[pathway.listOfSpecies(i).x pathway.listOfSpecies(i).y...
        pathway.listOfSpecies(i).w pathway.listOfSpecies(i).h];
    if strcmpi(pathway.listOfSpecies(i).type, 'SIMPLE_MOLECULE')
        metRec=[metRec;position];
        metText=[metText;pathway.listOfSpecies(i).name];
        metX=[metX;position(1)+0.3*position(3)];
        metY=[metY;position(2)+0.5*position(4)];
    end
    if strcmpi(pathway.listOfSpecies(i).type, 'PROTEIN')
        %The color of the enzyme can be specified in the pathway object.
        %The default color is used if there is no 'color' field present.
        %If neither of the fluxes are above the cutoff value, then use the
        %default color
        faceColor=defaultColor;
        if isfield(pathway.listOfSpecies(i),'color')
            if any(pathway.listOfSpecies(i).color)
                %Check that the fluxes are specified
                if isfield(pathway.listOfSpecies(i),'flux') && isfield(pathway.listOfSpecies(i),'referenceFlux')
                    %Check that the value of one of the fluxes is above the cutoff value
                    if abs(pathway.listOfSpecies(i).referenceFlux)>=cutOff || abs(pathway.listOfSpecies(i).flux)>=cutOff
                        faceColor=pathway.listOfSpecies(i).color;
                    end
                end
            end
        end
        
        %If the reaction is associated with a sign change and either of the
        %fluxes are larger than the cutoff value, then use a different frame
        edgeColor=[0 0 0];
        if isfield(pathway.listOfSpecies(i),'signChange')
            if any(pathway.listOfSpecies(i).signChange)
                if (pathway.listOfSpecies(i).signChange==true) && (abs(pathway.listOfSpecies(i).referenceFlux)>=cutOff || abs(pathway.listOfSpecies(i).flux)>=cutOff)
                    edgeColor=[1 0.6 0.2];
                end
            end
        end
        
        rxnFace=[rxnFace;faceColor];
        rxnEdge=[rxnEdge;edgeColor];
        
        %Draw smaller boxes if expression data should also be printed
        if drawExpression==true
            rxnRec=[rxnRec;position-[0 0 40 0]];
            
            %Draw rectangles representing each gene. Max 8 circles are used.
            %Their width is hardcoded and so is the width of the area where
            %they are printed (40)
            if any(pathway.listOfSpecies(i).expA)
                %Calculate the position of the dots. 4 in each column.
                firstFourX=pathway.listOfSpecies(i).x+pathway.listOfSpecies(i).w-34;
                fourY=pathway.listOfSpecies(i).y:pathway.listOfSpecies(i).h/4:pathway.listOfSpecies(i).y+pathway.listOfSpecies(i).h;
                
                nExp=numel(pathway.listOfSpecies(i).expA);
                
                %One log10-fold change is total up or down
                colorCodes=getColorCodes(pathway.listOfSpecies(i).expA, pathway.listOfSpecies(i).expB);
                
                %Plot the circles
                for j=1:min(nExp,4)
                    expRec=[expRec;firstFourX fourY(j) 16 pathway.listOfSpecies(i).h/4];
                    expFace=[expFace;colorCodes{j}];
                end
                
                %If a second row is needed
                if nExp>4
                    secondFourX=pathway.listOfSpecies(i).x+pathway.listOfSpecies(i).w-16;
                    for j=5:min(nExp,8)
                        expRec=[expRec;secondFourX fourY(j-4) 16 pathway.listOfSpecies(i).h/4];
                        expFace=[expFace;colorCodes{j}];
                    end
                end
            end
        else
            rxnRec=[rxnRec;position];
        end
        
        %If no fluxes are specified then should only the name be printed.
        %Flux values should only be printed if at least one of the fluxes
        %is above the cutoff value
        textToPrint={};
        
        %NOTE:  This is done since I know that I have names including '_' which is a reserved
        %       character in LaTeX. This should be done more systematically
        
        %NOTE:  The notation used means that only flux values with with less than
        %       nine characters can be used
        
        textToPrint{1,1}=strrep(pathway.listOfSpecies(i).name,'_','\_');
        
        %Check that the fluxes are specified
        %NOTE: Should you only print if there are two fluxes?
        if isfield(pathway.listOfSpecies(i),'flux') && isfield(pathway.listOfSpecies(i),'referenceFlux')
            %Check that the value of one of the fluxes is above the cutoff value
            if ~isempty(pathway.listOfSpecies(i).flux) && ~isempty(pathway.listOfSpecies(i).referenceFlux)
                if abs(pathway.listOfSpecies(i).referenceFlux)>=cutOff || abs(pathway.listOfSpecies(i).flux)>=cutOff
                    textToPrint{2,1}=num2str(pathway.listOfSpecies(i).flux,'%0.4g');
                    textToPrint{3,1}=num2str(pathway.listOfSpecies(i).referenceFlux,'%0.4g');
                end
            end
        end
        
        xPos=ones(numel(textToPrint),1)*(position(1)+3);
        distance=position(4)/(numel(textToPrint)+1);
        yPos=position(2)+distance:distance:(position(2)+position(4)-distance);
        
        rxnText=[rxnText;textToPrint];
        rxnX=[rxnX;xPos];
        rxnY=[rxnY;yPos'];
    end
end

%*** Print everything ***

%Compartment rectangles
for i=1:size(compRec,1)
    rectangle('edgecolor',[1 0.8 0.2], 'facecolor', [1 0.91 0.55], 'linewidth', 2, 'curvature', [0.05 0.05],...
        'position',compRec(i,:));
    rectangle('edgecolor',[1 0.8 0.2], 'facecolor', [1 1 1], 'linewidth', 1, 'curvature', [0.05 0.05],...
        'position',compRec(i,:)+[12 12 -24 -24]);
end

%Reaction lines
line(rxnLinesX',rxnLinesY','linewidth',0.1,'color',[0 0 0]);

%Arrow lines
hold on;
for i=1:size(rxnArrows,1)
    plotArrow(rxnArrows(i,1),rxnArrows(i,2),rxnArrows(i,3),rxnArrows(i,4),'linewidth',0.1,'edgecolor',rxnArrowsColor(i,:),'facecolor',rxnArrowsColor(i,:));
end
hold off;

%Enzyme lines
line(enzLinesX',enzLinesY','marker','o','markeredgecolor',[0 0 0],'markersize',0.5,...
    'markerfacecolor',[1 0.969 0.922],'linewidth',0.1,'color',[0 0 0]);

%Metabolite rectangles
for i=1:size(metRec,1)
    rectangle('edgecolor', [0 0 0], 'facecolor', [1 1 0.78], 'curvature', [0.8 0.8],...
        'position',metRec(i,:),'linewidth',0.1);
end

%Reaction rectangles
for i=1:size(rxnRec,1)
    rectangle('edgecolor', rxnEdge(i,:), 'facecolor', rxnFace(i,:), 'linewidth',0.1,...
        'curvature', [0.1 0.1], 'position',rxnRec(i,:));
end

%Expression rectangles
for i=1:size(expRec,1)
    rectangle('facecolor', expFace(i,:), 'linewidth',0.1, 'curvature', [0.1 0.1],...
        'position',expRec(i,:));
end

%Metabolite text
%Temporary thing to get it right for small networks as well. Should be
%calculated better!
if dimension(3)<5000
    textSize=3;
else
    textSize=0.5;
end

text(metX,metY,metText,'fontname','Small Fonts','fontsize',textSize,'interpreter','tex',...
    'HorizontalAlignment','center','verticalalignment','middle');

%Reaction text
text(rxnX,rxnY,rxnText,'fontname','Small Fonts','fontsize',textSize,'interpreter','tex',...
    'verticalalignment','middle','HorizontalAlignment','left');

%Compartment text. The text is centered and at the bottom of the compartment
text(compX, compY,compText,'fontname','Small Fonts','fontsize',7,'interpreter','none',...
    'HorizontalAlignment','center','verticalalignment','middle');
end

function handles = plotArrow( x1,y1,x2,y2,varargin )
%
% plotArrow - plots an arrow to the current plot
%
% format:   handles = plotArrow( x1,y1,x2,y2 [,options...] )
%
% input:    x1,y1   - starting point
%           x2,y2   - end point
%           options - come as pairs of "property","value" as defined for "line" and "patch"
%                     controls, see matlab help for listing of these properties.
%                     note that not all properties where added, one might add them at the end of this file.
%
%                     additional options are:
%                     'headwidth':  relative to complete arrow size, default value is 0.07
%                     'headheight': relative to complete arrow size, default value is 0.15
%                     (encoded are maximal values if pixels, for the case that the arrow is very long)
%
% output:   handles - handles of the graphical elements building the arrow
%
% Example:  plotArrow( -1,-1,15,12,'linewidth',2,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5] );
%           plotArrow( 0,0,5,4,'linewidth',2,'headwidth',0.25,'headheight',0.33 );
%           plotArrow;   % will launch demo

% =============================================
% constants (can be edited)
% =============================================
%NOTE: These should be 100 times smaller
alpha       = 15;   % head length
beta        = 7;   % head width
max_length  = 11;
%max_width   = 7;

% =============================================
% check if head properties are given
% =============================================
% if ratio is always fixed, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        switch lower(varargin{c*2-1})
            % head properties - do nothing, since handled above already
            case 'headheight',alpha = max( min( varargin{c*2},1 ),0.01 );
            case 'headwidth', beta = max( min( varargin{c*2},1 ),0.01 );
        end
    end
end

% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;                                % make sure no devision by zero occurs
teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
cs          = cos(teta);                                    % rotation matrix
ss          = sin(teta);
R           = [cs -ss;ss cs];
line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
head_length = min( line_length*alpha,max_length );
head_width  = min( line_length*beta,max_length );
x0          = x2*cs + y2*ss;                                % build head coordinats
y0          = -x2*ss + y2*cs;
coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
h1          = plot( [x1,x2],[y1,y2],'k' );
h2          = patch( coords(1,:),coords(2,:),[0 0 0] );

% =============================================
% return handles
% =============================================
handles = [h1 h2];

% =============================================
% check if styling is required
% =============================================
% if no styling, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        switch lower(varargin{c*2-1})
            
            % only patch properties
            case 'edgecolor',   set( h2,'EdgeColor',varargin{c*2} );
            case 'facecolor',   set( h2,'FaceColor',varargin{c*2} );
            case 'facelighting',set( h2,'FaceLighting',varargin{c*2} );
            case 'edgelighting',set( h2,'EdgeLighting',varargin{c*2} );
                
                % only line properties
            case 'color'    , set( h1,'Color',varargin{c*2} );
                
                % shared properties
            case 'linestyle', set( handles,'LineStyle',varargin{c*2} );
            case 'linewidth', set( handles,'LineWidth',varargin{c*2} );
            case 'parent',    set( handles,'parent',varargin{c*2} );
                
                % head properties - do nothing, since handled above already
            case 'headwidth',
            case 'headheight',
                
        end
    end
end
end
