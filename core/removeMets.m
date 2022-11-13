function reducedModel=removeMets(model,metsToRemove,isNames,removeUnusedRxns,removeUnusedGenes,removeUnusedComps)
% removeMets
%   Deletes a set of metabolites from a model
%
%   model             a model structure
%   metsToRemove      either a cell array of metabolite IDs, a logical vector
%                     with the same number of elements as metabolites in the model,
%                     of a vector of indexes to remove
%   isNames           true if the supplied mets represent metabolite names
%                     (as opposed to IDs). This is a way to delete
%                     metabolites in several compartments at once without
%                     knowing the exact IDs. This only works if metsToRemove
%                     is a cell array (opt, default false)
%   removeUnusedRxns  remove reactions that are no longer in use (opt,
%                     default false)
%   removeUnusedGenes remove genes that are no longer in use (opt,
%                     default false)
%   removeUnusedComps remove compartments that are no longer in use (opt,
%                     default false)
%
%   reducedModel      an updated model structure
%
%   Usage: reducedModel=removeMets(model,metsToRemove,isNames,...
%           removeUnusedRxns,removeUnusedGenes,removeUnusedComps)
if ~islogical(metsToRemove) && ~isnumeric(metsToRemove)
    metsToRemove=convertCharArray(metsToRemove);
end

if nargin<3
    isNames=false;
end

if nargin<4
    removeUnusedRxns=false;
end

if nargin<5
    removeUnusedGenes=false;
end

if nargin<6
    removeUnusedComps=false;
end

%Check that metsToRemove is a cell array
if isNames==true && ~iscell(metsToRemove)
    error('Must supply a cell array of strings if isNames=true');
end

reducedModel=model;

if isNames==false
    indexesToDelete=getIndexes(model,metsToRemove,'mets');
else
    indexesToDelete=[];
    for i=1:numel(metsToRemove)
        indexesToDelete=[indexesToDelete;find(strcmp(metsToRemove(i),model.metNames))];
    end
end

%Remove metabolites
if ~isempty(indexesToDelete)
    reducedModel.mets(indexesToDelete)=[];
    reducedModel.S(indexesToDelete,:)=[];
    if isfield(reducedModel,'b')
        reducedModel.b(indexesToDelete,:)=[];
    end
    if isfield(reducedModel,'metNames')
        reducedModel.metNames(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metComps')
        reducedModel.metComps(indexesToDelete)=[];
    end
    if isfield(reducedModel,'inchis')
        reducedModel.inchis(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metSmiles')
        reducedModel.metSmiles(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metFormulas')
        reducedModel.metFormulas(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metMiriams')
        reducedModel.metMiriams(indexesToDelete)=[];
    end
    if isfield(reducedModel,'unconstrained')
        reducedModel.unconstrained(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metFrom')
        reducedModel.metFrom(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metCharges')
        reducedModel.metCharges(indexesToDelete)=[];
    end
end

%Remove unused reactions
if removeUnusedRxns==true
    %Get unused reactions
    [~, a]=find(reducedModel.S);
    rxnsToRemove=1:numel(reducedModel.rxns);
    rxnsToRemove(a)=[];
    reducedModel=removeReactions(reducedModel,rxnsToRemove,false,removeUnusedGenes);
end

%Remove unused compartments
if removeUnusedComps==true
    oldComps=reducedModel.comps;
    I=ismember(1:numel(oldComps),reducedModel.metComps);
    if ~all(I)
        reducedModel.comps(~I)=[];
        reducedModel.compNames(~I)=[];
        if isfield(reducedModel,'compOutside')
            reducedModel.compOutside(~I)=[];
        end
        if isfield(reducedModel,'compMiriams')
            reducedModel.compMiriams(~I)=[];
        end
        [~, J]=ismember(oldComps(reducedModel.metComps),reducedModel.comps);
        reducedModel.metComps=J;
    end
end
end
