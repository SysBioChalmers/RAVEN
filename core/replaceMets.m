function [model, removedRxns, idxDuplRxns]=replaceMets(model,metabolite,replacement,verbose,identifiers)
% replaceMets
%   Replaces metabolite names and annotation with replacement metabolite
%   that is already in the model. If this results in duplicate metabolites,
%   the replacement metabolite will be kept, while the S matrix is updated
%   to use the replacement metabolite instead. At the end, contractModel is
%   run to remove any duplicate reactions that might have occured.
%
% Input:
%   model           a model structure
%   metabolite      string with name of metabolite to be replace
%   replacement     string with name of replacement metabolite
%   verbose         logical whether to print the ids of reactions that
%                   involve the replaced metabolite (optional, default
%                   false)
%   identifiers     true if 'metabolite' and 'replacement' refer to
%                   metabolite identifiers instead of metabolite names
%                   (optional, default false)
% 
% Output:
%   model           model structure with selected metabolites replaced
%   removedRxns     identifiers of duplicate reactions that were removed
%   idxDuplRxns     indices of removedRxns in the input model
%
% Note: This function is useful when the model contains both 'oxygen' and
% 'o2' as metabolite names. If 'oxygen' and 'o2' are identifiers instead,
% then the 'identifiers' flag should be set to true.
%
% Usage: [model, removedRxns, idxDuplRxns] = replaceMets(model, metabolite, replacement, verbose)

metabolite=char(metabolite);
replacement=char(replacement);

if nargin<4 || isempty(verbose)
    verbose=false;
end
if nargin<5
    identifiers = false;
end

% Find occurence of replacement metabolites. Annotation will be taken from
% first metabolite found.
if identifiers
    repIdx = find(strcmp(replacement,model.mets));
else
    repIdx = find(strcmp(replacement,model.metNames));
end
if isempty(repIdx)
    error('The replacement metabolite cannot be found in the model.');
end

% Change name and information from metabolite to replacement metabolite
if identifiers
    metIdx = find(strcmp(metabolite,model.mets));
else
    metIdx = find(strcmp(metabolite,model.metNames));
end
if isempty(metIdx)
    error('The to-be-replaced metabolite cannot be found in the model.');
end

[~, rxnsWithMet] = find(model.S(metIdx,:));
rxnsWithMet = unique(rxnsWithMet);
if verbose==true
    fprintf('\n\nThe following reactions involve the to-be-replaced metabolite:\n')
    fprintf('%s\n',model.rxns{rxnsWithMet})
    fprintf('\n')
end

model.metNames(metIdx) = model.metNames(repIdx(1));
if isfield(model,'metFormulas')
    model.metFormulas(metIdx) = model.metFormulas(repIdx(1));
end
if isfield(model,'metMiriams')
    model.metMiriams(metIdx) = model.metMiriams(repIdx(1));
end
if isfield(model,'metCharges')
    model.metCharges(metIdx) = model.metCharges(repIdx(1));
end
if isfield(model,'metDeltaG')
    model.metDeltaG(metIdx) = model.metDeltaG(repIdx(1));
end
if isfield(model,'inchis')
    model.inchis(metIdx) = model.inchis(repIdx(1));
end
if isfield(model,'metSmiles')
    model.metSmiles(metIdx) = model.metSmiles(repIdx(1));
end

idxDelete=[];
if identifiers
    if any(ismember(metIdx,repIdx))
        error('The to-be-replaced metabolite and the replacement metabolite are the same.');
    end
    originalStoch = model.S(metIdx,rxnsWithMet);
    model.S(repIdx,rxnsWithMet) = model.S(repIdx,rxnsWithMet) + originalStoch;
    model.S(metIdx,rxnsWithMet) = 0;
    idxDelete = metIdx;
else
    % Run through replacement metabolites and their compartments. If any of the
    % to-be-replaced metabolites is already present (checked by
    % metaboliteName[compartment], then the replacement metabolite is kept and
    % the to-be-replace metabolite ID deleted.
    
    % Build list of metaboliteName[compartment]
    comps = model.comps(:);
    metCompsN = strcat(lower(model.metNames(:)),'[',comps(model.metComps(:)),']');
    
    % Within each metaboliteName[compartment] group, the first replacement
    % metabolite is kept and the other members are merged into it.
    groups = unique(metCompsN(repIdx),'stable');
    for i = 1:numel(groups)
        metCompsNidx = find(strcmp(groups{i}, metCompsN));
        keepIdx = metCompsNidx(find(ismember(metCompsNidx,repIdx),1));
        others = metCompsNidx(metCompsNidx ~= keepIdx);
        if ~isempty(others)
            model.S(keepIdx,:) = model.S(keepIdx,:) + sum(model.S(others,:),1);
            idxDelete = [idxDelete; others]; % Make list of metabolite IDs to delete
        end
    end
end
if ~isempty(idxDelete)
    repMets = model.mets(repIdx);
    model = removeMets(model,idxDelete);
    % Indices of the replacement metabolites shift when metabolites with a
    % lower index are deleted.
    [~, repIdx] = ismember(repMets,model.mets);
    repIdx(repIdx==0) = [];
end

% This could now have created duplicate reactions. Contract model.
if isempty(idxDelete)
    removedRxns = cell(0,1);
    idxDuplRxns = zeros(0,1);
else
    originalRxns = model.rxns;
    [model, removedRxns] = contractModel(model,[],repIdx);
    [~, idxDuplRxns] = ismember(removedRxns,originalRxns);
end
end
