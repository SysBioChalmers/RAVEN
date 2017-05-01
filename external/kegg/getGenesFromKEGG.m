function model=getGenesFromKEGG(keggPath,koList)
% getGenesFromKEGG
%   Retrieves information on all genes stored in KEGG database
%
%   keggPath    if keggGenes.mat is not in the RAVEN\external\kegg directory,
%               this function will attempt to read data from a local FTP dump
%               of the KEGG database. keggPath is the path to the root of
%               this database
%
%   koList      the number of genes in KEGG is very large. koList can be a
%               cell array with KO identifiers, in which case only genes
%               belonging to one of those KEGG orthologies are retrieved (opt,
%               default all KOs with associated reactions)
%
%   model       a model structure generated from the database. The following
%               fields are filled
%               id:             'KEGG'
%               description:    'Automatically generated from KEGG database'
%               rxns:           KO ids
%               rxnNames:       Name for each entry
%               genes:          IDs for all the genes. Genes are saved as
%                               organism abbreviation:id (same as in KEGG).
%                               'HSA:124' for example is alcohol dehydrogenase
%                               in Homo sapiens
%               rxnGeneMat      A binary matrix that indicates whether a
%                               specific gene is present in a KO id
%
%   If the file keggGenes.mat is in the RAVEN\external\kegg directory it
%   will be loaded instead of parsing of the KEGG files. If it does not exist
%   it will be saved after parsing of the KEGG files. In general, you should
%   remove the keggGenes.mat file if you want to rebuild the model structure
%   from a newer version of KEGG.
%
%   Usage: model=getGenesFromKEGG(keggPath,koList)
%
%   Simonas Marcisauskas, 2017-05-02
%

%NOTE: This is how one entry looks in the file

% ENTRY       K11440                      KO
% NAME        gbsB
% DEFINITION  choline dehydrogenase [EC:1.1.1.1]
% CLASS       Metabolism; Amino Acid Metabolism; Glycine, serine and threonine
%             metabolism [PATH:ko00260]
% DBLINKS     RN: R08557 R08558
%             COG: COG1454
% GENES       BSU: BSU31050(gbsB)
%             BCY: Bcer98_1477
%             BLI: BL02563(gbsB)
%             BLD: BLi03269(gbsB)
%             BCL: ABC0979(gbsB)
%             BAY: RBAM_028150(gbsB)
%             BPU: BPUM_2734(gbsB)
% ADDITIONAL  INFO...
% ///

%The file is not tab-delimited. Instead each label is 12 characters
%(except for '///')

%Check if the genes have been parsed before and saved. If so, load the
%model.
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
genesFile=fullfile(ravenPath,'external','kegg','keggGenes.mat');
if exist(genesFile, 'file')
    fprintf(['NOTE: Importing KEGG genes from ' strrep(genesFile,'\','/') '.\n']);
    load(genesFile);
else
    fprintf(['Cannot locate ' strrep(genesFile,'\','/') ' and will try to generate it from the local KEGG database.\n']);
    if ~exist(fullfile(keggPath,'ko'),'file') || ~exist(fullfile(keggPath,'reaction'),'file')
        EM=fprintf(['The files ''ko'' and ''reaction'' cannot be located at ' strrep(keggPath,'\','/') '/ and should be downloaded from the KEGG FTP.\n']);
        dispEM(EM);
    else
        %Get all KOs that are associated to reactions
        allKOs=getAllKOs(keggPath);
   
        %Since the list of genes will be accessed many times I use a hash table
        geneMap=containers.Map('KeyType','char','ValueType','int32');
   
        %Add new functionality in the order specified in models
        model.id='KEGG';
        model.description='Automatically generated from KEGG database';
   
        %Preallocate memory
        model.rxns=cell(numel(allKOs),1);
        model.rxnNames=cell(numel(allKOs),1);
   
        %Reserve space for 500000 genes
        model.genes=cell(500000,1);
   
        %Load information on KO ID, name, and associated genes
        fid = fopen(fullfile(keggPath,'ko'), 'r');
   
        %Keeps track of how many KOs that have been added
        koCounter=0;
        nGenes=0;
        addingGenes=false;
   
        %These contain the mapping between KOs and genes and are used to
        %generate the model.rxnGeneMat in the end
        koIndex=zeros(5000000,1);
        geneIndex=koIndex;
        nMappings=0;
   
        skipRead=false;
   
        %Loop through the file
        while 1
          %Get the next line
          if skipRead==false
            tline = fgetl(fid);
          else
            skipRead=false;
          end
   
          %Abort at end of file
          if ~ischar(tline)
              break;
          end
   
          %Skip '///'
          if numel(tline)<12
              addingGenes=false;
              continue;
          end
   
          %Check if it's a new reaction
          if strcmp(tline(1:12),'ENTRY       ')
              %Check if it should be added
              koID=tline(13:18);
   
              if ismember(koID,allKOs)
                 addMe=true;
                 koCounter=koCounter+1;
              else
                 addMe=false;
                 continue;
              end
   
              %Add reaction ID (always 6 characters)
              model.rxns{koCounter}=koID;
              model.rxnNames{koCounter}=''; %Will be overwritten if it exists
          end
   
          %To get here we must be in a KO that should be added
          if addMe==true
              %Add name
              if strcmp(tline(1:12),'DEFINITION  ')
                  model.rxnNames{koCounter}=tline(13:end);
              end
   
              %Check if it's time to start adding genes
              if strcmp(tline(1:12),'GENES       ')
                  addingGenes=true;
              end
   
              %Add genes
              if addingGenes==true
                 %We are now adding genes for the current KO. All gene rows start
                 %with 12 spaces. If this is not true it means that there is
                 %additional annotation such as AUTHORS. This is not common but it
                 %exists.
                 if strcmp(tline(1:12),'            ') || strcmp(tline(1:12),'GENES       ')
                     geneLine=tline(13:end);
   
                     %Check if the line is from a new organism of from the same as
                     %before
                     if strcmp(geneLine(4),':')
                         % if organism id contains 3 letters;
                         currentOrganism=geneLine(1:3);
                         %Parse the string
                         genes=regexp(geneLine(6:end),' ','split');
                         genes=strcat(currentOrganism,':',genes(:));
                     else if strcmp(geneLine(5),':')
                         % if organism id contains 4 letters;
                         currentOrganism=geneLine(1:4);
                         %Parse the string
                         genes=regexp(geneLine(7:end),' ','split');
                         genes=strcat(currentOrganism,':',genes(:));
                         end
                     end
   
                     %Add the genes to the gene list
                     for i=1:numel(genes)
                         geneString=genes{i};
                         if geneMap.isKey(geneString);
                            index=geneMap(geneString);
                         else
                            nGenes=nGenes+1;
                            model.genes{nGenes}=geneString;
                            index=nGenes;
                            geneMap(geneString)=index;
                         end
   
                         %Add the connection to the KO
                         nMappings=nMappings+1;
                         koIndex(nMappings)=koCounter;
                         geneIndex(nMappings)=index;
                     end
                 else
                     %Now we want to throw away everything until the next
                     %entry
                     while 1
                        tline = fgetl(fid);
                        if strcmp(tline,'///')
                            %When the new entry is found, skip reading one line to
                            %fit with the rest of the code
                            skipRead=true;
                            break;
                        end
                     end
                 end
              end
          end
        end
        %Close the file
        fclose(fid);
   
        %If too much space was allocated, shrink the model
        model.rxns=model.rxns(1:koCounter);
        model.rxnNames=model.rxnNames(1:koCounter);
        model.genes=model.genes(1:nGenes);
   
        %Retrieve and add the gene associations
        model.rxnGeneMat=sparse(koIndex(1:nMappings),geneIndex(1:nMappings),ones(nMappings,1));
   
        %To make sure the size is correct if the last KOs don't have genes
        if size(model.rxnGeneMat,1)~=koCounter
            model.rxnGeneMat(koCounter,1)=0;
        end
   
        %Save the model structure
        save(genesFile,'model');
    end
end

%Only get the KOs in koList
I=~ismember(model.rxns,koList);
model=removeReactions(model,I,true,true);
end

function allKOs=getAllKOs(keggPath)
    %Retrieves all KOs that are associated to reactions. This is because
    %the number of genes in KEGG is very large so without this parsing it
    %would take many hours

    allKOs={};

    %First check if the reactions have already been parsed
    [ST, I]=dbstack('-completenames');
    ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
    rxnsFile=fullfile(ravenPath,'external','kegg','keggRxns.mat');
    if exist(rxnsFile, 'file')
        fprintf(['NOTE: Importing KEGG ORTHOLOGY list from ' strrep(rxnsFile,'\','/') '.\n']);
        load(rxnsFile,'model');

        %Loop through the reactions and add the corresponding genes
        for i=1:numel(model.rxns)
            if isstruct(model.rxnMiriams{i})
                %Get all KOs
                allKOs=[allKOs;model.rxnMiriams{i}.value(strcmpi(model.rxnMiriams{i}.name,'kegg.orthology'))];
            end
        end
    else
        %Parse the reaction file instead
        %First load information on reaction ID, reaction name, KO, pathway, and ec-number
        fid = fopen(fullfile(keggPath,'reaction'), 'r');
        orthology=false;
        while 1
          %Get the next line
          tline = fgetl(fid);

          %Abort at end of file
          if ~ischar(tline)
              break;
          end

          %Skip '///'
          if numel(tline)<12
              continue;
          end

          %Check if it's a new reaction
          if strcmp(tline(1:12),'ENTRY       ')
              orthology=false;
          end

          if strcmp(tline(1:9),'REFERENCE')
              orthology=false;
          end

          %Add KO-ids
          if numel(tline)>16
              if strcmp(tline(1:16),'ORTHOLOGY   KO: ') || strcmp(tline(1:16),'            KO: ') || strcmp(tline(1:12),'ORTHOLOGY   ') || orthology==true
                  if strcmp(tline(13:16),'KO:') %This is in the old version
                    allKOs=[allKOs;tline(17:22)];
                  else
                      if strcmp(tline(1:12),'ORTHOLOGY   ')
                        %This means that it found one KO in the new format and that
                        %subsequent lines might be other KOs
                        orthology=true;
                      end
                      allKOs=[allKOs;tline(13:18)];
                  end
              end
          end
        end

        %Close the file
        fclose(fid);
    end
    allKOs=unique(allKOs);
end
