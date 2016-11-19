function phylDistStruct=getPhylDist(keggPath,onlyInKingdom)
% getPhylDist
%   Calculates distance between species in KEGG based on systematic name.
%
%   keggPath        pathway to the KEGG database files
%   onlyInKingdom   if true, it generates a distance matrix with distance 
%                   Inf for organisms from another domains
%                   (Prokaryota, Eukaryota) (opt, default false)
%
%   phylDistStruct  a structure with a list of organism ids and a matrix
%                   that specifies their pairwise distances
%
%   This simple metric is based on the number of nodes two organisms are
%   away from each other in KEGG.
%
%   Usage: phylDistStruct=getPhylDist(keggPath,onlyInKingdom)
%
%   Rasmus Agren, 2012-12-17
%

if nargin<2
    onlyInKingdom=false;
end

%Check if the reactions have been parsed before and saved. If so, load the
%model.
[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));
% Adding escape characters, if some parent folders contain spaces or
% exclamation marks (for Unix systems). For Windows, all the parent folders
% are just put between the double quotation brackets
if isunix
    ravenPath = regexprep(ravenPath,'\ ','\\ ');
    ravenPath = regexprep(ravenPath,'\!','\\!');
elseif ispc
    for i=1:(length(strfind(ravenPath,'\')))
        if i==1
            ravenPath = regexprep(ravenPath,'\\','\\"',i);
        elseif i==length(strfind(ravenPath,'\'))
            ravenPath = regexprep(ravenPath,'\\','"\\',i);    
        else
            ravenPath = regexprep(ravenPath,'\\','"\\"',i);
        end
    end
end

distFile=fullfile(ravenPath,'external','kegg','keggPhylDist.mat');
if exist(distFile, 'file')
    fprintf(['NOTE: Importing KEGG phylogenetic distance matrix from ' strrep(distFile,'\','/') '.\n']);
    load(distFile);
else
    %Download required files from KEGG if it doesn't exist in the directory
    downloadKEGG(keggPath);

    %Open the file that describes the naming of the species
    fid = fopen(fullfile(keggPath,'taxonomy'), 'r');

    phylDistStruct.ids={};

    %Keeps the categories for each organism
    orgCat={};

    currentCat={}; %Keeps track of the current category
    depth=0; %Keeps track of the current tree depth

    %Loop through the file
    orgCounter=0;
    while 1
      %Get the next line
      tline = fgetl(fid);

      %Abort at end of file
      if ~ischar(tline)
          break;
      end

      if any(tline)
          %Check if it's a new category
          if tline(1)=='#'
              %Find the first space (=depth +1)
              sPos=strfind(tline,' ')-1; %Should always exist
              sPos=sPos(1);

              %If we have stepped back one step in the tree
              if sPos<depth
                  currentCat=currentCat(1:sPos);
              end
              depth=sPos;

              currentCat{depth}=tline(sPos+2:end);
          else
              orgCounter=orgCounter+1;
              %It is an organism
              %Get the id between first and second white space
              sPos=find(isstrprop(tline, 'wspace')); %Should always exist
              phylDistStruct.ids{orgCounter}=tline(sPos(1)+1:sPos(2)-1);
              orgCat{orgCounter}=currentCat;
          end
      end
    end

    %Generate a distance matrix (very straight forward here, not neat)
    phylDistStruct.distMat=zeros(numel(phylDistStruct.ids));
    for i=1:numel(phylDistStruct.ids)
        for j=1:numel(phylDistStruct.ids)
            if onlyInKingdom==true
               if ~strcmp(orgCat{i}(1),orgCat{j}(1))
                   phylDistStruct.distMat(i,j)=Inf;
                   continue;
               end
            end
            %Calculate the distance between then
            dist=numel(orgCat{i})-numel(orgCat{j});
            if dist>0
                aCat=orgCat{i}(1:end-dist);
            else
                aCat=orgCat{i};
            end
            if dist<0
                bCat=orgCat{j}(1:end+dist);
            else
                bCat=orgCat{j};
            end

            %Loop through the categories and stop when they are the same
            for k=numel(aCat):-1:1
               if strcmp(aCat{k},bCat{k})
                  break; 
               end
            end
            phylDistStruct.distMat(i,j)=dist+numel(aCat)-k;
        end
    end
    %Save the structure
    save(distFile,'phylDistStruct');
end
end