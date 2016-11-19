function pathway = constructPathwayFromCelldesigner(inputFile)
% constructPathwayFromCelldesigner
%	Constructs a pathway structure from a CellDesigner XML file.
%
%	inputFile   string representing the pathway of the Celldesigner XML 
%               file
%
%	pathway
%       listOfCompartments      structure with information on the compartments
%                               in the Celldesigner file
%           x                   the x position
%           y                   the y position
%           h                   the height
%           w                   the width
%           compartment         string representing the compartment
%           id                  alias of the compartment. Used if there are
%                               several instances of the same compartment
%           name                the name of the compartment
%       listOfSpecies           structure with information on each species in
%                               the Celldesigner file. Both metabolites and
%                               enzymes are species
%           x                   the x position
%           y                   the y position
%           h                   the height
%           w                   the width
%           species             string representing the species
%           alias               alias of the species. Used if the species is
%                               present in several places in the map
%           name                the name of the species
%           type                the type of the species. At the moment only
%                               'SIMPLE_MOLECULE' (for metabolites) and
%                               'PROTEIN' (for enzymes) are supported
%           note                string that can be used on enzymes to link them
%                               to the corresponding reaction in a model
%       listOfReactions         structure with information on each of the
%                               reactions in the Celldesigner file
%           reversible          true if the reaction is reversible
%           middlePoint         vector with the x and y points of the middle
%                               point of the reaction. All reactions are
%                               defines as having one reactant and one product.
%                               They may also have a number of modifiers (such
%                               as cofactors or enzymes). Each of these
%                               modifiers are drawn as being connected to the
%                               middle point. The middle point is calculated
%                               from the positions of the reactant and the
%                               product
%           componentList       structure with information on the reactant,
%                               product, and modifiers of the reaction
%               alias           string representing the alias of the
%                               species. Corresponds to the alias in
%                               listOfSpecies
%               species         string representing the species.
%                               Corresponds to the species in listOfSpecies
%               anchor          vector with the x and y points of the
%                               connecting point of the species
%               baseReaction    true if the species is the base reactant
%               baseProduct     true if the species is the base product
%               toArrow         true if the line connecting the species
%                               to the middle point should start with an
%                               arrow
%               fromArrow       true if the line connecting the species
%                               to the middle point should end with an
%                               arrow
%               type            'METABOLITE' or 'ENZYME'
%
%	Usage: pathway = constructPathwayFromCelldesigner(inputFile)
%
%   Rasmus Agren, 2010-12-16
%

%Loads the specified xml file using XML Toolbox
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

%Current path and xml toolbox path
cp=pwd;
xt=fullfile(ravenPath,'software','xml_toolbox');

%Change path, run the script and change back
cd(xt);
v = xml_parseany(fileread(inputFile));
cd(cp);

%Saves information on each compartment
for i=1:length(v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}.compartmentAlias)
   pathway.listOfCompartments(i).x=str2double(v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.bounds{1,1}.ATTRIBUTE.x);
   pathway.listOfCompartments(i).y=str2double(v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.bounds{1,1}.ATTRIBUTE.y);
   pathway.listOfCompartments(i).h=str2double(v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.bounds{1,1}.ATTRIBUTE.h);
   pathway.listOfCompartments(i).w=str2double(v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.bounds{1,1}.ATTRIBUTE.w);
   pathway.listOfCompartments(i).compartment=v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.ATTRIBUTE.compartment;
   pathway.listOfCompartments(i).id=v.model{1,1}.annotation{1,1}.listOfCompartmentAliases{1,1}...
       .compartmentAlias{1,i}.ATTRIBUTE.id;
   
  %Finds the name of the compartment
  for j=1:length(v.model{1,1}.listOfCompartments{1,1}.compartment)
     if strcmpi(v.model{1,1}.listOfCompartments{1,1}.compartment{1,j}.ATTRIBUTE.id,...
             pathway.listOfCompartments(i).compartment)
        pathway.listOfCompartments(i).name=v.model{1,1}.listOfCompartments{1,1}.compartment{1,j}.ATTRIBUTE.name;
        break;
     end
  end
end

%Saves information on each species and enzyme
for i=1:length(v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}.speciesAlias)
   pathway.listOfSpecies(i).x=str2double(v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.bounds{1,1}.ATTRIBUTE.x);
   pathway.listOfSpecies(i).y=str2double(v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.bounds{1,1}.ATTRIBUTE.y);
   pathway.listOfSpecies(i).h=str2double(v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.bounds{1,1}.ATTRIBUTE.h);
   pathway.listOfSpecies(i).w=str2double(v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.bounds{1,1}.ATTRIBUTE.w);
   pathway.listOfSpecies(i).alias=v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.ATTRIBUTE.id;  
   pathway.listOfSpecies(i).species=v.model{1,1}.annotation{1,1}.listOfSpeciesAliases{1,1}...
       .speciesAlias{1,i}.ATTRIBUTE.species;
  %Find the name and type of the species/enzyme
  for j=1:length(v.model{1,1}.listOfSpecies{1,1}.species)
     if strcmpi(v.model{1,1}.listOfSpecies{1,1}.species{1,j}.ATTRIBUTE.id,...
             pathway.listOfSpecies(i).species)
        pathway.listOfSpecies(i).name=v.model{1,1}.listOfSpecies{1,1}.species{1,j}.ATTRIBUTE.name;
        pathway.listOfSpecies(i).type=v.model{1,1}.listOfSpecies{1,1}.species{1,j}.annotation{1,1}...
            .speciesIdentity{1,1}.class{1,1}.CONTENT;
        
        %The Celldesigner chart can be linked to a model (e.g. for
        %visualization of fluxes) by making a note in Celldesigner for the
        %protein in interest (this can only be done for proteins at the
        %moment). The note (not protein note!) should only contain the
        %reaction id in the model which corresponds to the protein 
        %in Celldesigner.
        if strcmpi(pathway.listOfSpecies(i).type,'PROTEIN')
           if isfield(v.model{1,1}.listOfSpecies{1,1}.species{1,j},'notes')
               pathway.listOfSpecies(i).note=cellstr(strtrim(v.model{1,1}.listOfSpecies{1,1}.species{1,j}...
                   .notes{1,1}.html{1,1}.body{1,1}.CONTENT));
           end
        end
        break;
     end
  end
end

%Saves information on each reaction
for i=1:length(v.model{1,1}.listOfReactions{1,1}.reaction)
    %This is because it is standard for a reaction to be defined as
    %reversible
    if isfield(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.ATTRIBUTE,'reversible')
        pathway.listOfReactions(i).reversible=v.model{1,1}.listOfReactions{1,1}...
            .reaction{1,i}.ATTRIBUTE.reversible;
    else
        pathway.listOfReactions(i).reversible='true';
    end
    
    %NOTE: As far as I know there can only be one base reactant/product per
    %reaction. This method is written with that assumption.
    
    %Finds the alias of the base reactant
    baseReactant=v.model{1,1}.listOfReactions{1,1}...
        .reaction{1,i}.annotation{1,1}.baseReactants{1,1}.baseReactant{1,1}...
        .ATTRIBUTE.alias;
    
    %Finds the alias of the base product
    baseProduct=v.model{1,1}.listOfReactions{1,1}...
        .reaction{1,i}.annotation{1,1}.baseProducts{1,1}.baseProduct{1,1}...
        .ATTRIBUTE.alias;
    
    %Saves the middle point
    [x1,y1]=getBindingPos(pathway,baseReactant,...
        v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
        .baseReactants{1,1}.baseReactant{1,1}.linkAnchor{1,1}.ATTRIBUTE.position);
    [x2,y2]=getBindingPos(pathway,baseProduct,...
        v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
        .baseProducts{1,1}.baseProduct{1,1}.linkAnchor{1,1}.ATTRIBUTE.position);
    pathway.listOfReactions(i).middlePoint(1)=x1+(x2-x1)/2;
    pathway.listOfReactions(i).middlePoint(2)=y1+(y2-y1)/2;
    
    %Saves information on each of the components
    counter=1; %Keeps track of where in componentList to add a component
    
    %Adds the base reactant
    pathway.listOfReactions(i).componentList(counter).alias=v.model{1,1}...
        .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.baseReactants{1,1}...
        .baseReactant{1,1}.ATTRIBUTE.alias;
    pathway.listOfReactions(i).componentList(counter).species=v.model{1,1}...
        .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.baseReactants{1,1}...
        .baseReactant{1,1}.ATTRIBUTE.species;
    [x,y]=getBindingPos(pathway,baseReactant,v.model{1,1}.listOfReactions{1,1}...
        .reaction{1,i}.annotation{1,1}.baseReactants{1,1}.baseReactant{1,1}...
        .linkAnchor{1,1}.ATTRIBUTE.position);
    pathway.listOfReactions(i).componentList(counter).anchor(1)=x;
    pathway.listOfReactions(i).componentList(counter).anchor(2)=y;
    pathway.listOfReactions(i).componentList(counter).baseReactant='true';
    pathway.listOfReactions(i).componentList(counter).baseProduct='false';
    
    %toArrow should be true if the reaction is reversible
    pathway.listOfReactions(i).componentList(counter).toArrow=...
        pathway.listOfReactions(i).reversible;
    pathway.listOfReactions(i).componentList(counter).fromArrow='false';
    pathway.listOfReactions(i).componentList(counter).type='METABOLITE';
    counter=counter+1;
    
    %Adds the base product
    pathway.listOfReactions(i).componentList(counter).alias=v.model{1,1}...
        .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.baseProducts{1,1}...
        .baseProduct{1,1}.ATTRIBUTE.alias;
    pathway.listOfReactions(i).componentList(counter).species=v.model{1,1}...
        .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.baseProducts{1,1}...
        .baseProduct{1,1}.ATTRIBUTE.species;
    [x,y]=getBindingPos(pathway,baseProduct,v.model{1,1}.listOfReactions{1,1}...
        .reaction{1,i}.annotation{1,1}.baseProducts{1,1}.baseProduct{1,1}...
        .linkAnchor{1,1}.ATTRIBUTE.position);
    pathway.listOfReactions(i).componentList(counter).anchor(1)=x;
    pathway.listOfReactions(i).componentList(counter).anchor(2)=y;    
    pathway.listOfReactions(i).componentList(counter).baseReactant='false';
    pathway.listOfReactions(i).componentList(counter).baseProduct='true';
    pathway.listOfReactions(i).componentList(counter).toArrow='true';
    pathway.listOfReactions(i).componentList(counter).fromArrow='false';
    pathway.listOfReactions(i).componentList(counter).type='METABOLITE';
    counter=counter+1;
        
    %Adds the non-base reactants
    %Not all reactions have any non-base reactants
    if isfield(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
            ,'listOfReactantLinks')
        for j=1:length(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfReactantLinks{1,1}.reactantLink)
            pathway.listOfReactions(i).componentList(counter).alias=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfReactantLinks{1,1}.reactantLink{1,j}.ATTRIBUTE.alias;
            pathway.listOfReactions(i).componentList(counter).species=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfReactantLinks{1,1}.reactantLink{1,j}.ATTRIBUTE.reactant;
            [x,y]=getBindingPos(pathway,pathway.listOfReactions(i).componentList(counter).alias...
                ,v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfReactantLinks{1,1}.reactantLink{1,j}.linkAnchor{1,1}...
                .ATTRIBUTE.position);
            pathway.listOfReactions(i).componentList(counter).anchor(1)=x;
            pathway.listOfReactions(i).componentList(counter).anchor(2)=y;
            pathway.listOfReactions(i).componentList(counter).baseReactant='false';
            pathway.listOfReactions(i).componentList(counter).baseProduct='false';
            pathway.listOfReactions(i).componentList(counter).toArrow='false';
            pathway.listOfReactions(i).componentList(counter).fromArrow='false';
            pathway.listOfReactions(i).componentList(counter).type='METABOLITE';
            counter=counter+1;
        end
    end
    
    %Adds the non-base products
    %Not all reactions have any non-base products
    if isfield(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
            ,'listOfProductLinks')
        for j=1:length(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfProductLinks{1,1}.productLink)
            pathway.listOfReactions(i).componentList(counter).alias=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfProductLinks{1,1}.productLink{1,j}.ATTRIBUTE.alias;
            pathway.listOfReactions(i).componentList(counter).species=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfProductLinks{1,1}.productLink{1,j}.ATTRIBUTE.product;
            [x,y]=getBindingPos(pathway,pathway.listOfReactions(i).componentList(counter).alias...
                ,v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfProductLinks{1,1}.productLink{1,j}.linkAnchor{1,1}...
                .ATTRIBUTE.position);
            pathway.listOfReactions(i).componentList(counter).anchor(1)=x;
            pathway.listOfReactions(i).componentList(counter).anchor(2)=y;
            pathway.listOfReactions(i).componentList(counter).baseReactant='false';
            pathway.listOfReactions(i).componentList(counter).baseProduct='false';
            pathway.listOfReactions(i).componentList(counter).toArrow='true';
            pathway.listOfReactions(i).componentList(counter).fromArrow='false';
            pathway.listOfReactions(i).componentList(counter).type='METABOLITE';
            counter=counter+1;
        end
    end
    
    %Adds the modifiers (enzymes)
    %Not all reactions have any modifiers
    if isfield(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
            ,'listOfModification')
        for j=1:length(v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfModification{1,1}.modification)
            pathway.listOfReactions(i).componentList(counter).alias=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.listOfModification{1,1}...
                .modification{1,j}.linkTarget{1,1}.ATTRIBUTE.alias;
            pathway.listOfReactions(i).componentList(counter).species=v.model{1,1}...
                .listOfReactions{1,1}.reaction{1,i}.annotation{1,1}.listOfModification{1,1}...
                .modification{1,j}.linkTarget{1,1}.ATTRIBUTE.species;
            [x,y]=getBindingPos(pathway,pathway.listOfReactions(i).componentList(counter).alias...
                ,v.model{1,1}.listOfReactions{1,1}.reaction{1,i}.annotation{1,1}...
                .listOfModification{1,1}.modification{1,j}.linkTarget{1,1}...
                .linkAnchor{1,1}.ATTRIBUTE.position);
            pathway.listOfReactions(i).componentList(counter).anchor(1)=x;
            pathway.listOfReactions(i).componentList(counter).anchor(2)=y;
            pathway.listOfReactions(i).componentList(counter).baseReactant='false';
            pathway.listOfReactions(i).componentList(counter).baseProduct='false';
            pathway.listOfReactions(i).componentList(counter).toArrow='false';
            pathway.listOfReactions(i).componentList(counter).fromArrow='false';
            pathway.listOfReactions(i).componentList(counter).type='ENZYME';
            counter=counter+1;
        end
    end
end

function [x,y]=getBindingPos(pathway,speciesAlias,bindingSite)
% getBindingPos
%   Calculates the position of a given binding site for a specified species
%
%   pathway         the pathway structure which should be modified
%   speciesAlias    the alias of the species
%   bindingSite     The binding site in Celldesigner is a string of up to
%                   three letters that defines which of the 16 bindings 
%                   sites on each species to use. They are written as 
%                   N (north), NNW (north north west) and so on
%
%   [x, y]          The position of the binding site for the species
%
%   Usage: [x,y]=getBindingPos(pathway,speciesAlias,bindingSite)
%
%   Rasmus Ågren, 2010-12-16
%

%Find the species
for i=1:length(pathway.listOfSpecies)
   if strcmpi(pathway.listOfSpecies(i).alias,speciesAlias)
       xpos=pathway.listOfSpecies(i).x;
       ypos=pathway.listOfSpecies(i).y;
       h=pathway.listOfSpecies(i).h;
       w=pathway.listOfSpecies(i).w;
       
       %Species marked as "PROTEIN" should be drawn as rectangles. There is
       %one binding site in each corner and three in between on each side.
       if strcmpi(pathway.listOfSpecies(i).type,'PROTEIN')
           %The binding sites clockwise from the top position
           if strcmpi(bindingSite,'N')
              x=xpos+0.5*w;
              y=ypos;
              break;
           end
           if strcmpi(bindingSite,'NNE')
              x=xpos+0.75*w;
              y=ypos;
              break;
           end
           if strcmpi(bindingSite,'NE')
              x=xpos+w;
              y=ypos;
              break;
           end
           if strcmpi(bindingSite,'ENE')
              x=xpos+w;
              y=ypos+0.25*h;
              break;
           end
           if strcmpi(bindingSite,'E')
              x=xpos+w;
              y=ypos+0.5*h;
              break;
           end
           if strcmpi(bindingSite,'ESE')
              x=xpos+w;
              y=ypos+0.75*h;
              break;
           end
           if strcmpi(bindingSite,'SE')
              x=xpos+w;
              y=ypos+h;
              break;
           end
           if strcmpi(bindingSite,'SSE')
              x=xpos+0.75*w;
              y=ypos+h;
              break;
           end
           if strcmpi(bindingSite,'S')
              x=xpos+0.5*w;
              y=ypos+h;
              break;
           end
           if strcmpi(bindingSite,'SSW')
              x=xpos+0.25*w;
              y=ypos+h;
              break;
           end
           if strcmpi(bindingSite,'SW')
              x=xpos;
              y=ypos+h;
              break;
           end
           if strcmpi(bindingSite,'WSW')
              x=xpos;
              y=ypos+0.75*h;
              break;
           end
           if strcmpi(bindingSite,'W')
              x=xpos;
              y=ypos+0.5*h;
              break;
           end
           if strcmpi(bindingSite,'WNW')
              x=xpos;
              y=ypos+0.25*h;
              break;
           end
           if strcmpi(bindingSite,'NW')
              x=xpos;
              y=ypos;
              break;
           end
           if strcmpi(bindingSite,'NNW')
              x=xpos+0.25*w;
              y=ypos;
              break;
           end
       end
       
       %Species marked as "SIMPLE_MOLECULE" should be drawn as ellipses. 
       %There is one binding site in each extreme point and three in 
       %between on each side
       if strcmpi(pathway.listOfSpecies(i).type,'SIMPLE_MOLECULE')
           %Find the center of the ellipse
           centerX=xpos+0.5*w;
           centerY=ypos+0.5*h;
           
           %The binding sites clockwise from the top position
           if strcmpi(bindingSite,'N')
              angle=(4/8)*pi;
           end
           if strcmpi(bindingSite,'NNE')
              angle=(5/8)*pi;
           end
           if strcmpi(bindingSite,'NE')
              angle=(6/8)*pi;
           end
           if strcmpi(bindingSite,'ENE')
              angle=(7/8)*pi;
           end
           if strcmpi(bindingSite,'E')
              angle=(8/8)*pi;
           end
           if strcmpi(bindingSite,'ESE')
              angle=(9/8)*pi;
           end
           if strcmpi(bindingSite,'SE')
              angle=(10/8)*pi;
           end
           if strcmpi(bindingSite,'SSE')
              angle=(11/8)*pi;
           end
           if strcmpi(bindingSite,'S')
              angle=(12/8)*pi;
           end
           if strcmpi(bindingSite,'SSW')
              angle=(13/8)*pi;
           end
           if strcmpi(bindingSite,'SW')
              angle=(14/8)*pi;
           end
           if strcmpi(bindingSite,'WSW')
              angle=(15/8)*pi;
           end
           if strcmpi(bindingSite,'W')
              angle=(16/8)*pi;
           end
           if strcmpi(bindingSite,'WNW')
              angle=(1/8)*pi;
           end
           if strcmpi(bindingSite,'NW')
              angle=(2/8)*pi;
           end
           if strcmpi(bindingSite,'NNW')
              angle=(3/8)*pi;
           end
           x=centerX-(0.5*w)*cos(angle);
           y=centerY-(0.5*h)*sin(angle);
           break;
       end
   end
end