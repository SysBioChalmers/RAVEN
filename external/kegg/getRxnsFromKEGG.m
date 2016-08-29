function model=getRxnsFromKEGG(keggPath,keepUndefinedStoich,keepIncomplete, keepGeneral)
% getRxnsFromKEGG
%   Retrieves information on all reactions stored in KEGG database
%
%   keggPath            this function reads data from a local FTP dump of 
%                       the KEGG database. keggPath is the pathway to the
%                       root of the database
%   keepUndefinedStoich include reactions in the form n A <=> n+1 A. These
%                       will be dealt with as two separate metabolites 
%                       (opt, default true)
%   keepIncomplete      include reactions which have been labelled as
%                       "incomplete", "erroneous" or "unclear" (opt,
%                       default true)
%   keepGeneral         include reactions which have been labelled as
%                       "general reaction". These are reactions on the form
%                       "an aldehyde <=> an alcohol", and are therefore
%                       unsuited for modelling purposes. Note that not all
%                       reactions have this type of annotation, and the
%                       script will therefore not be able to remove all
%                       such reactions (opt, default false)
%
%   model     a model structure generated from the database. The following
%             fields are filled
%             id:             'KEGG'
%             description:    'Automatically generated from KEGG database'
%             rxns:           KEGG reaction ids
%             rxnNames:       Name for each reaction entry
%             mets:           KEGG compound ids. If the equations use
%                             stoichiometry such as ID(n+1) then the whole 
%                             expression is saved as the id
%             eccodes:        Corresponding ec-number if available
%             rxnMiriams:     Contains reaction specific information such as
%                             KO id and pathways that the reaction is 
%                             associated to
%             S:              Stoichiometric matrix
%             lb:             -1000 for all reactions
%             ub:             1000 for all reactions
%             rev:            1 for reversible and 0 for irreversible. For
%                             reactions present in pathway maps the reversibility 
%                             is taken from there
%             b:              0 for all metabolites
%
%   Reactions on the form A <=> A + B will not be loaded. If the file
%   keggRxns.mat is in the RAVEN directory it will be loaded instead of
%   parsing of the KEGG files. If it does not exist it will be saved after
%   parsing of the KEGG files. In general, you should remove the
%   keggRxns.mat file if you want to rebuild the model structure from a
%   newer version of KEGG.
%
%   Usage: model=getRxnsFromKEGG(keggPath,keepUndefinedStoich,keepIncomplete,keepGeneral)
%
%   Rasmus Agren, 2012-12-16
%

%NOTE: This is how one entry looks in the file

% ENTRY       R00001                      Reaction
% NAME        Polyphosphate polyphosphohydrolase
% DEFINITION  Polyphosphate + n H2O <=> (n+1) Oligophosphate
% EQUATION    C00890 + n C00001 <=> (n+1) C02174
% ENZYME      3.6.1.10
% ///

%The file is not tab-delimited. Instead each label is 12 characters
%(except for '///')

if nargin<2
    keepUndefinedStoich=true;
end
if nargin<3
    keepIncomplete=true;
end
if nargin<4
    keepGeneral=true;
end

%Check if the reactions have been parsed before and saved. If so, load the
%model.
[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
rxnsFile=fullfile(ravenPath,'external','kegg','keggRxns.mat');
if exist(rxnsFile, 'file')
    fprintf(['NOTE: Importing KEGG reactions from ' strrep(rxnsFile,'\','/') '.\n']);
    load(rxnsFile);
else
    %Download required files from KEGG if it doesn't exist in the directory
    downloadKEGG(keggPath);
    
    %Add new functionality in the order specified in models
    model.id='KEGG';
    model.description='Automatically generated from KEGG database';

    %Preallocate memory for 10000 reactions
    model.rxns=cell(10000,1);
    model.rxnNames=cell(10000,1);
    model.eccodes=cell(10000,1);
    model.subSystems=cell(10000,1);
    model.rxnMiriams=cell(10000,1);
    equations=cell(10000,1); %Temporarily stores the equations
    isIncomplete=false(10000,1);
    isGeneral=false(10000,1);

    %First load information on reaction ID, reaction name, KO, pathway, and ec-number
    fid = fopen(fullfile(keggPath,'reaction'), 'r');

    %Keeps track of how many reactions that have been added
    rxnCounter=0;

    %Loop through the file
    orthology=false;
    pathway=false;
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
          rxnCounter=rxnCounter+1;

          %Add empty strings where there should be such
          model.rxnNames{rxnCounter}='';
          model.eccodes{rxnCounter}='';
          model.subSystems{rxnCounter}='';
          equations{rxnCounter}='';

          %Add reaction ID (always 6 characters)
          model.rxns{rxnCounter}=tline(13:18);
          orthology=false;
          pathway=false;
      end

      %Add name
      if strcmp(tline(1:12),'NAME        ')
          model.rxnNames{rxnCounter}=tline(13:end);
      end
      
      %Add whether the comment includes "incomplete", "erroneous" or "unclear"
      if strcmp(tline(1:12),'COMMENT     ')
          %Read all text until '///' or 'RPAIR'
          commentText=tline(13:end);
          while 1
            tline = fgetl(fid);
            if ~strcmp(tline(1:3),'///') && ~strcmp(tline(1:3),'RPA') && ~strcmp(tline(1:3),'ENZ')
               commentText=strcat(commentText,' ',tline);
            else
                break;
            end
          end
          upperLine=upper(commentText);
          if any(strfind(upperLine,'INCOMPLETE')) || any(strfind(upperLine,'ERRONEOUS')) || any(strfind(upperLine,'UNCLEAR'))
            isIncomplete(rxnCounter)=true;
          end
          if any(strfind(upperLine,'GENERAL REACTION')==1) %It should start this way
            isGeneral(rxnCounter)=true;
          end
              
          %Go to next iteration if it is '///'
          if numel(tline)<12
              continue;
          end
      end
      
      %Add ec-number
      if strcmp(tline(1:12),'ENZYME      ')
          model.eccodes{rxnCounter}=tline(13:end);
      end
      if numel(tline)>8
        if strcmp(tline(1:9),'REFERENCE')
            pathway=false;
            orthology=false;
        end
      end

      %Add KO-ids
      if numel(tline)>16
          if strcmp(tline(1:16),'ORTHOLOGY   KO: ') || strcmp(tline(1:16),'            KO: ') || strcmp(tline(1:12),'ORTHOLOGY   ') || orthology==true
              pathway=false;
              %Check if KO has been added already (each reaction may belong to
              %several)
              if isstruct(model.rxnMiriams{rxnCounter})
                  addToIndex=numel(model.rxnMiriams{rxnCounter}.name)+1;
              else
                  addToIndex=1;
              end

              tempStruct=model.rxnMiriams{rxnCounter};
              tempStruct.name{addToIndex,1}='urn:miriam:kegg.ko'; %WARNING: This is not a real MIRIAM identifier
              if strcmp(tline(13:16),'KO:') %This is in the old version
                tempStruct.value{addToIndex,1}=tline(17:22);
              else
                  if strcmp(tline(1:12),'ORTHOLOGY   ')
                    %This means that it found one KO in the new format and that
                    %subsequent lines might be other KOs
                    orthology=true;
                  end
                  tempStruct.value{addToIndex,1}=tline(13:18);  
              end
              model.rxnMiriams{rxnCounter}=tempStruct;
          end
      end

      %Add pathways
      if numel(tline)>18
          if strcmp(tline(1:18),'PATHWAY     PATH: ') || strcmp(tline(1:18),'            PATH: ') || strcmp(tline(1:12),'PATHWAY     ') || pathway==true
              %Check if annotation has been added already
              if isstruct(model.rxnMiriams{rxnCounter})
                  addToIndex=numel(model.rxnMiriams{rxnCounter}.name)+1;
              else
                  addToIndex=1;
              end

              tempStruct=model.rxnMiriams{rxnCounter};
              tempStruct.name{addToIndex,1}='urn:miriam:kegg.pathway';
              %If it's the old version
              if strcmp(tline(14:17),'PATH:')
                tempStruct.value{addToIndex,1}=tline(19:25);
              else
                %If it's the new version
                tempStruct.value{addToIndex,1}=tline(13:19);
                pathway=true;
              end

              %Don't do this if the pathway is rn01100 (Metabolic pathways)
              if ~strcmp('rn01100',tempStruct.value{addToIndex,1})
                model.rxnMiriams{rxnCounter}=tempStruct;

                %Also save the subSystems entry as being the first path found
                if ~any(model.subSystems{rxnCounter})
                    if strcmp(tline(14:17),'PATH:')
                        model.subSystems{rxnCounter}=tline(28:end);
                    else
                        model.subSystems{rxnCounter}=tline(22:end);
                    end
                end
              end
          end
      end
    end

    %Close the file
    fclose(fid);

    %This is done here since the the indexes won't match since some reactions
    %are removed along the way
    isIncomplete=model.rxns(isIncomplete);
    isGeneral=model.rxns(isGeneral);

    %If too much space was allocated, shrink the model
    model.rxns=model.rxns(1:rxnCounter);
    model.rxnNames=model.rxnNames(1:rxnCounter);
    model.eccodes=model.eccodes(1:rxnCounter);
    equations=equations(1:rxnCounter);
    model.rxnMiriams=model.rxnMiriams(1:rxnCounter);
    model.subSystems=model.subSystems(1:rxnCounter);

    %Then load the equations from another file. This is because the equations
    %are easier to retrieve from there

    %The format is rxnID: equation
    %The reactions should have been loaded in the exact same order
    fid = fopen(fullfile(keggPath,'reaction.lst'), 'r');

    %Loop through the file
    for i=1:rxnCounter
      %Get the next line
      tline = fgetl(fid);

      equations{i}=tline(9:end);
    end

    %Close the file
    fclose(fid);

    %Construct the S matrix and list of metabolites
    [S mets badRxns]=constructS(equations);
    model.S=S;
    model.mets=mets;

    %There is some limited evidence for directionality in
    %reaction_mapformula.lst. The information there concerns how the reactions
    %are drawn in the KEGG maps. If a reaction is irreversible in the same
    %direction for all maps, then I consider is irreversible, otherwise
    %reversible. Also, not all reactions are present in the maps, so not all
    %will have directionality. They will be considered to be reversible.

    %The format is R00005: 00330: C01010 => C00011
    %Generate a reversibility structure with the fields:
    %rxns: reaction ids
    %product: one met id that is a product. This is because the reactions
    %might be written in another direction compared to in the reactions.lst
    %file
    %rev: 1 if reversible, otherwise 0
    reversibility.rxns={};
    reversibility.product={};
    reversibility.rev=[];

    fid = fopen(fullfile(keggPath,'reaction_mapformula.lst'), 'r');
    while 1
      %Get the next line
      tline = fgetl(fid);

      %Abort at end of file
      if ~ischar(tline)
          break;
      end

      rxn=tline(1:6);
      prod=tline(end-5:end);
      rev=any(strfind(tline,'<=>'));
      if isempty(reversibility.rxns)
        reversibility.rxns{1}=rxn;
        reversibility.product{1}=prod;
        reversibility.rev(1)=rev;
      else
        %Check if the reaction was added before. It's an ordered list, so only
        %check the last element
        if strcmp(reversibility.rxns(end),rxn)
            %If it's reversible in the new reaction or reversible in the old reaction
            %then set (keep) to be reversible
            if rev==1 || reversibility.rev(end)==1
                reversibility.rev(end)=1;
            else
                %This means that the reaction was already loaded, that it was
                %irreversible before and irreversible in the new reaction.
                %However, it could be that they are written in diferent
                %directions. If the product differ, then set to be reversible.
                %This assumes that the reactions are written with the same
                %metabolite as the last one if they are in the same direction.
                if ~strcmp(prod,reversibility.product(end))
                    reversibility.rev(end)=1;
                end
            end
        else
            reversibility.rxns=[reversibility.rxns;rxn];
            reversibility.product=[reversibility.product;prod];
            reversibility.rev=[reversibility.rev;rev];
        end
      end
    end
    fclose(fid);

    %Update the reversibility
    model.rev=ones(rxnCounter,1);
    %Match the reaction ids
    irrevIDs=find(reversibility.rev==0);
    [crap I]=ismember(reversibility.rxns(irrevIDs),model.rxns);
    [crap prodMetIDs]=ismember(reversibility.product(irrevIDs),model.mets);
    model.rev(I)=0;

    %See if the reactions are written in the same order in model.S
    linearInd=sub2ind(size(model.S), prodMetIDs, I);
    changeOrder=I(model.S(linearInd)<0);
    model.S(:,changeOrder)=model.S(:,changeOrder).*-1; %Change the order of these reactions

    %Add some stuff to get a correct model structure
    model.ub=ones(rxnCounter,1)*1000;
    model.lb=model.rev*-1000;
    model.c=zeros(rxnCounter,1);
    model.b=zeros(numel(model.mets),1);
    model=removeReactions(model,badRxns,true,true);
    
    %Save the model structure
    save(rxnsFile,'model','isGeneral','isIncomplete');
end

%Delete reaction which are labeled as "incomplete", "erroneous", "unclear"
%or "general reaction" (depending on settings.
if keepGeneral==false
    model=removeReactions(model,intersect(isGeneral,model.rxns),true,true);
end
if keepIncomplete==false
    model=removeReactions(model,intersect(isIncomplete,model.rxns),true,true);
end

%Delete reactions involving undefined stoichiometry. These metabolites have
%an ID containing the letter "n" or "m"
if keepUndefinedStoich==false
    I=cellfun(@any,strfind(model.mets,'n')) | cellfun(@any,strfind(model.mets,'m'));
    [crap J]=find(model.S(I,:));
    model=removeReactions(model,J,true,true);
end
end
