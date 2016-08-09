function newModel=addGenes(model,genesToAdd)
% addGenes
%   Adds genes to a model
%
%   model        a model structure
%   genesToAdd   the genes genesToAdd can have the following fields:
%                genes          cell array with unique strings that 
%                               identifies each gene. Only character which are
%                               allowed in SBML ids are allowed (mainly a-z,
%                               0-9 and '_'). However, there is no check
%                               for this performed, as it only matters if
%                               the model should be exported to SBML
%                geneShortNames cell array of gene abbreviations (opt,
%                               default '')
%                geneMiriams    cell array with MIRIAM structures (opt,
%                               default [])
%
%   newModel     an updated model structure
%
%   NOTE: This function does not make extensive checks about MIRIAM formats,
%   forbidden characters or such.
%
%   Usage: newModel=addGenes(model,genesToAdd)
%
%   Rasmus Agren, 2013-08-01
%

newModel=model;

if isempty(genesToAdd)
    return;
end

%Check some stuff regarding the required fields
if ~isfield(genesToAdd,'genes')
    dispEM('genes is a required field in genesToAdd');
end

if ~iscellstr(genesToAdd.genes)
    dispEM('genesToAdd.genes must be a cell array of strings');
end

illegalCells=regexp(genesToAdd.genes,'[();:]', 'once');
dispEM('Illegal character(s) in gene names:',true,genesToAdd.genes(~cellfun(@isempty,illegalCells)));

%Number of genes
nGenes=numel(genesToAdd.genes);
nOldGenes=numel(model.genes);
filler=cell(nGenes,1);
filler(:)={''};
largeFiller=cell(nOldGenes,1);
largeFiller(:)={''};

%Check that no gene ids are already present in the model
I=ismember(genesToAdd.genes,model.genes);
if any(I)
	dispEM('One or more elements in genesToAdd.genes are already present in model.genes');
else
    newModel.genes=[newModel.genes;genesToAdd.genes(:)];
end

%Some more checks and if they pass then add each field to the structure
if isfield(genesToAdd,'geneShortNames')
   if numel(genesToAdd.geneShortNames)~=nGenes
       dispEM('genesToAdd.geneShortNames must have the same number of elements as genesToAdd.genes');
   end
   if ~iscellstr(genesToAdd.geneShortNames)
        dispEM('genesToAdd.geneShortNames must be a cell array of strings');
   end
   %Add empty field if it doesn't exist
   if ~isfield(newModel,'geneShortNames')
        newModel.geneShortNames=largeFiller;
   end
   newModel.geneShortNames=[newModel.geneShortNames;genesToAdd.geneShortNames(:)];
else
    %Add empty strings if structure is in model
    if isfield(newModel,'geneShortNames')
       newModel.geneShortNames=[newModel.geneShortNames;filler];
    end
end

%Don't check the type of geneMiriams
if isfield(genesToAdd,'geneMiriams')
   if numel(genesToAdd.geneMiriams)~=nGenes
       dispEM('genesToAdd.geneMiriams must have the same number of elements as genesToAdd.genes');
   end
   %Add empty field if it doesn't exist
   if ~isfield(newModel,'geneMiriams')
        newModel.geneMiriams=cell(nOldGenes,1);
   end
   newModel.geneMiriams=[newModel.geneMiriams;genesToAdd.geneMiriams(:)]; 
else
    if isfield(newModel,'geneMiriams')
       newModel.geneMiriams=[newModel.geneMiriams;cell(nGenes,1)]; 
    end
end

if isfield(genesToAdd,'geneComps')
   if numel(genesToAdd.geneComps)~=nGenes
       dispEM('genesToAdd.geneComps must have the same number of elements as genesToAdd.genes');
   end
   %Add empty field if it doesn't exist
   if ~isfield(newModel,'geneComps')
        newModel.geneComps=ones(nOldGenes,1);
        dispEM('Adding genes with compartment information to a model without such information. All existing genes will be assigned to the first compartment',false);
   end
   newModel.geneComps=[newModel.geneComps;genesToAdd.geneComps(:)]; 
else
    if isfield(newModel,'geneComps')
       newModel.geneComps=[newModel.geneComps;ones(nGenes,1)]; 
       fprintf('NOTE: The added genes will be assigned to the first compartment\n');
    end
end

if isfield(newModel,'geneFrom')
    newModel.geneFrom=[newModel.geneFrom;filler];
end

if isfield(newModel,'rxnGeneMat')
    newModel.rxnGeneMat=[newModel.rxnGeneMat sparse(numel(model.rxns),nGenes)];
end
end
