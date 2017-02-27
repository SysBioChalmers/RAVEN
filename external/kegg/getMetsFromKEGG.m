function model=getMetsFromKEGG(keggPath)
% getMetsFromKEGG
%   Retrieves information on all metabolites stored in KEGG database
%
%   keggPath    if keggMets.mat is not in the RAVEN\external\kegg directory,
%               this function will attempt to read data from a local FTP dump
%               of the KEGG database. keggPath is the path to the root of
%               this database
%
%   model       a model structure generated from the database. The following
%               fields are filled
%               id:             'KEGG'
%               description:    'Automatically generated from KEGG database'
%               mets:           KEGG compound ids
%               metNames:       Compound name. Only the first name will be
%                               saved if there are several synonyms
%               metMiriams:     If there is a CHEBI id available, then that
%                               will be saved here
%               inchis:         InChI string for the metabolite
%               metFormulas:    The chemical composition of the metabolite.
%                               This will only be loaded if there is no InChI 
%                               string
%   If the file keggMets.mat is in the RAVEN\external\kegg directory it will be loaded 
%   instead of parsing of the KEGG files. If it does not exist it will be 
%   saved after parsing of the KEGG files. In general, you should remove the
%   keggMets.mat file if you want to rebuild the model structure from a
%   newer version of KEGG.
%               
%   Usage: model=getMetsFromKEGG(keggPath)
%
%   Eduard Kerkhoven, 2017-02-27
%

%NOTE: This is how one entry looks in the file

% ENTRY       C00001                      Compound
% NAME        H2O;
%             Water
% FORMULA     H2O
% MASS        18.0106
% REMARK      Same as: D00001
% REACTION    R00001 R00002 R00004 R00005 R00009 R00010 R00011 R00017
%             R00022 R00024 R00026 R00028 R00036 R00041 R00044 R00045
% ENZYME      1.1.1.160
% DBLINKS     PubChem: 7435
%             ChEBI: 29110

%Then a lot of info about the positions of the atoms and so on. It is not
%certain that each metabolite follows this structure exactly

%The file is not tab-delimited. Instead each label is 12 characters
%(except for '///')

%Check if the reactions have been parsed before and saved. If so, load the
%model.
[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
metsFile=fullfile(ravenPath,'external','kegg','keggMets.mat');
if exist(metsFile, 'file')
    fprintf(['NOTE: Importing KEGG metabolites from ' strrep(metsFile,'\','/') '.\n']);
    load(metsFile);
else
    fprintf(['Cannot locate ' strrep(metsFile,'\','/') '\nand will try to generate it from the local KEGG database.\n']);

    %Download required files from KEGG if it doesn't exist in the directory
    downloadKEGG(keggPath);
    
    %Add new functionality in the order specified in models
    model.id='KEGG';
    model.description='Automatically generated from KEGG database';

    %Preallocate memory for 20000 metabolites
    model.mets=cell(20000,1);
    model.metNames=cell(20000,1);
    model.metFormulas=cell(20000,1);
    model.metMiriams=cell(20000,1);

    %First load information on metabolite ID, metabolite name, composition, and
    %CHEBI
    fid = fopen(fullfile(keggPath,'compound'), 'r');

    %Keeps track of how many metabolites that have been added
    metCounter=0;

    %Loop through the file
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
          metCounter=metCounter+1;

          %Add empty strings where there should be such
          model.metNames{metCounter}='';
          model.metFormulas{metCounter}='';

          %Add compound ID (always 6 characters)
          model.mets{metCounter}=tline(13:18);
      end

      %Add name
      if strcmp(tline(1:12),'NAME        ')
          %If there are synonyms, then the last character is ';'
          if strcmp(tline(end),';')
                model.metNames{metCounter}=tline(13:end-1);
          else
                model.metNames{metCounter}=tline(13:end);
          end
      end

      %Add composition
      if strcmp(tline(1:12),'FORMULA     ')
          model.metFormulas{metCounter}=tline(13:end);
      end

      %Add CHEBI id
      if numel(tline)>19
          if strcmp(tline(1:19),'            ChEBI: ')
              chebiID=tline(20:end); %This is because there is sometimes more then one CHEBI index

              %Only load one id for now
              s=strfind(chebiID,' ');
              if any(s)
                 chebiID=chebiID(1:s-1);
              end
              miriamStruct.name{1}='obo.chebi:CHEBI';
              miriamStruct.value{1}=chebiID;
              model.metMiriams{metCounter}=miriamStruct;
          end
      end
    end

    %Close the file
    fclose(fid);

    %If too much space was allocated, shrink the model
    model.mets=model.mets(1:metCounter);
    model.metNames=model.metNames(1:metCounter);
    model.metFormulas=model.metFormulas(1:metCounter);
    model.metMiriams=model.metMiriams(1:metCounter);

    %If there was no CHEBI found, add the KEGG id as a metMiriams
    for i=1:numel(model.mets)
        if ~isstruct(model.metMiriams{i})
            miriamStruct.name{1}='kegg.compound';
            miriamStruct.value{1}=model.mets{i};
            model.metMiriams{i}=miriamStruct;
        end
    end

    %Then load the InChI strings from another file. Not all metabolites will be
    %present in the list

    inchIDs=cell(numel(model.mets),1);
    inchis=cell(numel(model.mets),1);

    %The format is metID*tab*string

    fid = fopen(fullfile(keggPath,'compound.inchi'), 'r');

    %Loop through the file
    counter=1;
    while 1
      %Get the next line
      tline = fgetl(fid);

      %Abort at end of file
      if ~ischar(tline)
          break;
      end

      %Get the ID and the InChI
      inchIDs{counter}=tline(1:6);
      inchis{counter}=tline(14:end);
      counter=counter+1;
    end

    %Close the file
    fclose(fid);

    inchIDs=inchIDs(1:counter-1);
    inchis=inchis(1:counter-1);

    %Find the metabolites that had InChI strings and add them to the model
    [a b]=ismember(inchIDs,model.mets);

    %If there were mets with InChIs but that were not in the list
    if ~all(a)
       dispEM('Not all metabolites with InChI strings were found in the original list');  
    end

    model.inchis=cell(numel(model.mets),1);
    model.inchis(:)={''};
    model.inchis(b)=inchis;

    %Remove composition if InChI was found
    model.metFormulas(b)={''};
    
    % Ensuring that all model.metMiriams.value consist only of strings, no
    % double
    for i=1:(numel(model.mets))
        for j=1:(numel(model.metMiriams{i}))
            if isa(model.metMiriams{i}.value{j},'double')
                model.metMiriams{i}.value{j}=num2str(model.metMiriams{i}.value{j});
            end
        end
    end
    
    %Saves the model
    save(metsFile,'model');
end
end
