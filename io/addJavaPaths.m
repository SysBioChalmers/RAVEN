% addJavaPaths
%   Adds the Apache POI classes to the static Java paths
%
%   Usage: addJavaPaths()
%
%   Rasmus Agren, 2016-02-19
%
function addJavaPaths()
    %Get the path to Apache POI
    [ST, I]=dbstack('-completenames');
    ravenPath=fileparts(fileparts(ST(I).file));
    poiPATH=fullfile(ravenPath,'software','apache-poi');

    toAdd={fullfile(poiPATH,'dom4j-1.6.1.jar');
        fullfile(poiPATH,'poi-3.8-20120326.jar');
        fullfile(poiPATH,'poi-ooxml-3.8-20120326.jar');
        fullfile(poiPATH,'poi-ooxml-schemas-3.8-20120326.jar');
        fullfile(poiPATH,'xmlbeans-2.3.0.jar');
        fullfile(poiPATH,'stax-api-1.0.1.jar')};

    %Open the javaclasspath.txt file or create it
    %otherwise   

    %Get the ones to add    
    for i=1:numel(toAdd)
        javaaddpath(toAdd{i});
    end
    %Throw an error to say that Matlab has to be restarted
end
