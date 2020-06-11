% addJavaPaths
%   Adds the Apache POI classes to the static Java paths
%
%   Usage: addJavaPaths()

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

existingPaths=javaclasspath();

%Add the paths that are not already present
for i=1:numel(toAdd)
    if ~any(ismember(existingPaths,toAdd{i}))
        javaaddpath(toAdd{i});
    end
end
end