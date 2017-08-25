function [colorCodes, signChange, errorFlag]= getColorCodes(referenceFluxes, fluxes, maxChange, defaultColor, upColor, downColor)
% getColorCodes
%	Calculates the coloring for a number of fluxes by comparing them to
%   reference fluxes.
%
%   referenceFluxes     vector of reference fluxes
%   fluxes              vector of fluxes. The number of elements in fluxes
%                       and referenceFluxes must be equal
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
%                       smaller than the reference flux (opt, default [1 0 0])

%   colorCodes          array list of colors in Matlab format in the same
%                       order as the fluxes
%   signChange          array list of boolean values where true indicates
%                       that there has been a sign change between the
%                       fluxes. Reactions with sign changes are not
%                       colored, but rather marked in another way. The
%                       order correspondes to the order of the fluxes
%   errorFlag           true if there has been an error
%
%   Usage: [colorCodes, signChange, errorFlag]=getColorCodes(referenceFluxes,...
%           fluxes, maxChange, defaultColor, upColor, downColor)
%
%   Rasmus Agren, 2010-12-16
%

if nargin<6
    downColor=[1 0 0];
end
if nargin<5
    upColor=[0 1 0];
end
if nargin<4
    defaultColor=[1 1 1];
end
if nargin<3
    maxChange=1;
end

%Checks that the flux vector and the reference flux vector have the same
%number of elements
if length(fluxes)~=length(referenceFluxes)
   colorCodes={};
   signChange={};
   errorFlag=1;
   fprintf('fluxes and referenceFluxes must have the same dimensions');
   return;
end

%Loops through the fluxes. If the reference flux is 0, then mark as
%upColor, if the flux is 0 then mark as downColor, and if both fluxes are 0
%then mark as defaultColor
for i=1:length(fluxes)
    signChange{i}=false;

    if referenceFluxes(i)==0 || fluxes(i)==0
        if referenceFluxes(i)==0 && fluxes(i)==0
            colorCodes{i}=defaultColor;
        else
            if referenceFluxes(i)==0
                colorCodes{i}=upColor;
            else
                colorCodes{i}=downColor;
            end
        end
    else
        %At the moment a negative flux that gets more negative is counted
        %as being upregulated. This might be counter intuitive.
        logvalue=log10(abs(fluxes(i))/abs(referenceFluxes(i)));

        %If there is a sign change
        if fluxes(i)*referenceFluxes(i)<0
            colorCodes{i}=defaultColor;
            signChange{i}=true;
        else
            colorCodes{i}=getColor(logvalue, defaultColor, upColor, downColor, maxChange);
        end
    end
end

function colorValue=getColor(logvalue, defaultColor, upColor, downColor, maxChange)
% getColor
%   Calculates the color for a specified logvalue.
%
%   logvalue            the logfold increase or desrease of a flux compared to the
%                       corresponding reference flux. Must be a positive
%                       value
%   defaultColor        a color in Matlab format to be used if there are no
%                       changes between the fluxes. This color is also used to
%                       calculate the transition between the colors for up and
%                       down regulated fluxes
%   upColor             a color in Matlab format to be used if the flux is
%                       larger than the reference flux
%   downColor           a color in Matlab format to be used if the flux is
%                       smaller than the reference flux
%   maxChange           the logfold increase or decrease that corresponds
%                       to full negative or full positive coloring
%
%   colorValue          vector with the calculated color
%
%   Usage: colorValue=getColor(logvalue, defaultColor, upColor, downColor, maxChange)
%
%   Rasmus Ågren, 2010-12-16
%

%If the flux has decreased
if logvalue<0
   %If the flux is lower than 10^-maxChange of the original then color as
   %downColor
   if logvalue<-maxChange
        colorValue=downColor;
   else
       %The color is linear from defaultColor to downColor
       colorValue=[defaultColor(1)+(downColor(1)-defaultColor(1))*logvalue/(-maxChange)...
                   defaultColor(2)+(downColor(2)-defaultColor(2))*logvalue/(-maxChange)...
                   defaultColor(3)+(downColor(3)-defaultColor(3))*logvalue/(-maxChange)];
   end
%If it has increased
else
   %If the flux is higher than 10^maxChange times the original then color green
   if logvalue>maxChange
        colorValue=upColor;
   else
       %The color is linear from defaultColor to upColor
       colorValue=[defaultColor(1)+(upColor(1)-defaultColor(1))*logvalue/(maxChange)...
                   defaultColor(2)+(upColor(2)-defaultColor(2))*logvalue/(maxChange)...
                   defaultColor(3)+(upColor(3)-defaultColor(3))*logvalue/(maxChange)];
   end
end
end