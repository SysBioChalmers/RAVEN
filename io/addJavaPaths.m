% addJavaPaths
%   Adds the Apache POI classes to the static Java paths
%
%   Usage: addJavaPaths()
%
%   Rasmus Agren, 2015-05-13
%
function addJavaPaths()
    %Get the path to Apache POI
    [ST, I]=dbstack('-completenames');
    ravenPath=fileparts(ST(I).file);
    poiPATH=fullfile(ravenPath,'software','apache-poi');
    
    toAdd={fullfile(poiPATH,'dom4j-1.6.1.jar');
        fullfile(poiPATH,'poi-3.8-20120326.jar');
        fullfile(poiPATH,'poi-ooxml-3.8-20120326.jar');
        fullfile(poiPATH,'poi-ooxml-schemas-3.8-20120326.jar');
        fullfile(poiPATH,'xmlbeans-2.3.0.jar');
        fullfile(poiPATH,'stax-api-1.0.1.jar')};
    
    %Open the javaclasspath.txt file or create it
    %otherwise
    fid=fopen(fullfile(prefdir,'javaclasspath.txt'),'r');
    if fid>0
        current=textscan(fid,'%s','Delimiter','\n');
        current=current{1};
        fclose(fid);
    else
        current={};
    end
    
    %Get the ones to add
    [~,I]=setdiff(upper(toAdd),upper(current));
    if any(I)
        fid=fopen(fullfile(prefdir,'javaclasspath.txt'),'at');
        fprintf(fid,'\n');
        for i=1:numel(I)
            fprintf(fid,[strrep(toAdd{I(i)},'\','\\') '\n']);
        end
        fclose(fid);
        
        %Throw an error to say that Matlab has to be restarted
        EM='RAVEN has added the Apache POI classes to the static Java path. Please restart Matlab to have these changes take effect';
        dispEM(EM);
    end 
end