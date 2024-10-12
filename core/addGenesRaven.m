function newModel=addGenesRaven(model,genesToAdd)
% addGenesRaven
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
%                geneShortNames cell array of gene abbreviations (optional,
%                               default '')
%                geneMiriams    cell array with MIRIAM structures (optional,
%                               default [])
%                proteins   cell array of protein names associated to
%                               each gene (optional, default '')
%
%   newModel     an updated model structure
%
%   NOTE: This function does not make extensive checks about MIRIAM formats,
%   forbidden characters or such.
%
% Usage: newModel=addGenesRaven(model,genesToAdd)

newModel=model;

if isempty(genesToAdd)
    return;
end

%Check some stuff regarding the required fields
if ~isfield(genesToAdd,'genes')
    EM='genes is a required field in genesToAdd';
    dispEM(EM);
else
    genesToAdd.genes=convertCharArray(genesToAdd.genes);
end

%Number of genes
nGenes=numel(genesToAdd.genes);
nOldGenes=numel(model.genes);
filler=cell(nGenes,1);
filler(:)={''};
largeFiller=cell(nOldGenes,1);
largeFiller(:)={''};

%Check that no gene ids are already present in the model
I=ismember(genesToAdd.genes,model.genes);
if all(I)
    warning('All genes in genesToAdd.genes are already present in model.genes');
    return
elseif any(I)
    existingGenes=strjoin(genesToAdd.genes(I), ', ');
    warning(['The following genes are already present in model.genes and will therefore not be added: ', existingGenes])
    genesToAdd.genes(I)=[];
    if isfield(genesToAdd,'geneShortNames')
        genesToAdd.geneShortNames(I)=[];
    end
    if isfield(genesToAdd,'proteins')
        genesToAdd.proteins(I)=[];
    end
    if isfield(genesToAdd,'geneMiriams')
        genesToAdd.geneMiriams(I)=[];
    end
else
    newModel.genes=[newModel.genes;genesToAdd.genes(:)];
end

%Some more checks and if they pass then add each field to the structure
if isfield(genesToAdd,'geneShortNames')
    genesToAdd.geneShortNames=convertCharArray(genesToAdd.geneShortNames);
    if numel(genesToAdd.geneShortNames)~=nGenes
        EM='genesToAdd.geneShortNames must have the same number of elements as genesToAdd.genes';
        dispEM(EM);
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
if isfield(genesToAdd,'proteins')
    genesToAdd.proteins=convertCharArray(genesToAdd.proteins);
    if numel(genesToAdd.proteins)~=nGenes
        EM='genesToAdd.proteins must have the same number of elements as genesToAdd.genes';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'proteins')
        newModel.proteins=largeFiller;
    end
    newModel.proteins=[newModel.proteins;genesToAdd.proteins(:)];
else
    %Add empty strings if structure is in model
    if isfield(newModel,'proteins')
        newModel.proteins=[newModel.proteins;filler];
    end
end


%Don't check the type of geneMiriams
if isfield(genesToAdd,'geneMiriams')
    if numel(genesToAdd.geneMiriams)~=nGenes
        EM='genesToAdd.geneMiriams must have the same number of elements as genesToAdd.genes';
        dispEM(EM);
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
        EM='genesToAdd.geneComps must have the same number of elements as genesToAdd.genes';
        dispEM(EM);
    end
    %Add empty field if it doesn't exist
    if ~isfield(newModel,'geneComps')
        newModel.geneComps=ones(nOldGenes,1);
        EM='Adding genes with compartment information to a model without such information. All existing genes will be assigned to the first compartment';
        dispEM(EM,false);
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
    newModel.rxnGeneMat=[newModel.rxnGeneMat,zeros(length(newModel.rxns),nGenes)];
end
end
