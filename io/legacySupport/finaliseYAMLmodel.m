function model = finaliseYAMLmodel(model, equations, enzStoich, isGECKO, verbose, ecGecko)
% finaliseYAMLmodel
%   Post-processing shared by both readYAMLmodel parsers (the
%   cobrapy-style canonical path in readYAMLmodel.m and the older
%   `!!omap` path in private/parseYAMLLegacy.m).
%
%   Fills the S-matrix from the per-rxn stoichiometry pairs, coerces
%   the cell-of-string numeric fields to numeric arrays, supplies
%   defaults for missing per-rxn / per-met / per-gene fields, drops
%   fields that ended up entirely empty, derives rxnGeneMat from
%   grRules, and runs the GECKO-specific finalisation when isGECKO
%   is true.
%
% Input:
%   model       the partially-built RAVEN model struct (lists and
%               cell-of-string fields populated by the parser).
%   equations   N x 3 cell array of (rxnIdx, metId, coeff) triples
%               that defines the stoichiometry matrix.
%   enzStoich   N x 3 cell array of (ecRxnIdx, enzymeId, coeff)
%               triples for ec.rxnEnzMat (only used when isGECKO).
%   isGECKO     true to apply ec.* finalisation.
%   verbose     true to print progress messages.
%   ecGecko     N x 2 cell array of (ecRxnIdx, ec-code) pairs from
%               the legacy multi-line eccodes path. Optional; not
%               used by the canonical parser.

if nargin < 6
    ecGecko = {};
end

if verbose
    fprintf('\nimporting completed\nfollow-up processing...');
end
[~, model.metComps] = ismember(model.metComps, model.comps);
[~, model.geneComps] = ismember(model.geneComps, model.comps);
[~, model.rxnComps] = ismember(model.rxnComps, model.comps);

% Fill S-matrix
if ~isempty(equations)
    rxnIdx = cellfun('isempty', equations(:,1));
    equations(rxnIdx,:) = '';
    if ~isempty(equations)
        rxnIdx = cell2mat(equations(:,1));
        [~,metIdx] = ismember(equations(:,2),model.mets);
        coeffs = cell2mat(equations(:,3));
        model.S=sparse(max(metIdx),max(rxnIdx));
        linearIndices = sub2ind([max(metIdx), max(rxnIdx)],metIdx,rxnIdx);
        model.S(linearIndices) = coeffs;
    end
end

% Convert strings to numeric
model.metCharges = str2double(model.metCharges);
model.lb = str2double(model.lb);
model.ub = str2double(model.ub);
model.rxnConfidenceScores = str2double(model.rxnConfidenceScores);
model.b = zeros(length(model.mets),1);
model.metDeltaG = str2double(model.metDeltaG);
model.rxnDeltaG = str2double(model.rxnDeltaG);

% Fill some other fields
model.annotation.defaultLB = min(model.lb);
model.annotation.defaultUB = max(model.ub);
if numel(model.lb)<numel(model.rxns) %No LB reported = min
    model.lb(end+1:numel(model.rxns)-numel(model.lb),1) = double(model.annotation.defaultLB);
end
if numel(model.ub)<numel(model.rxns) %No UB reported = max
    model.ub(end+1:numel(model.rxns)-numel(model.ub),1) = double(model.annotation.defaultUB);
end
if ~all(cellfun('isempty',model.rev))
    model.rev = str2double(model.rev);
else
    model.rev = [];
end
if numel(model.rev)<numel(model.rxns) %No rev reported, assume from LB and UB
    model.rev(end+1:numel(model.rxns)-numel(model.rev),1) = double(model.lb<0 & model.ub>0);
end

% Remove empty fields, otherwise fill to correct length
% Reactions
for i={'rxnNames','grRules','eccodes','rxnNotes','rxnReferences',...
       'rxnFrom','subSystems','rxnMiriams'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'rxns');
end
for i={'c'} % Zeros
   model = emptyOrFill(model,i{1},0,'rxns',true);
end
for i={'rxnConfidenceScores','rxnDeltaG'} % NaNs
   model = emptyOrFill(model,i{1},NaN,'rxns');
end
for i={'rxnComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'rxns');
end
% Metabolites
for i={'metNames','inchis','metFormulas','metMiriams','metFrom','metSmiles','metNotes'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'mets');
end
for i={'metCharges','unconstrained'} % Zeros
   model = emptyOrFill(model,i{1},0,'mets');
end
for i={'metDeltaG'} % NaNs
    model = emptyOrFill(model,i{1},NaN,'mets');
 end
for i={'metComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'mets');
end
% Genes
for i={'geneMiriams','geneShortNames','proteins'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'genes');
end
for i={'geneComps'} % Ones, assume first compartment
   model = emptyOrFill(model,i{1},1,'genes');
end
% Comps
for i={'compNames'} % Empty strings
   model = emptyOrFill(model,i{1},{''},'comps');
end
for i={'compOutside'} % First comp
   model = emptyOrFill(model,i{1},model.comps{1},'comps');
end

if isfield(model,'subSystems')
% If all entries are 1x1, then flatten
    if all(cellfun(@(x) numel(x) <= 1, model.subSystems))
        model.subSystems = transpose([model.subSystems{:}]);
    end
end

% Make rxnGeneMat fields and map to the existing model.genes field
if isfield(model,'grRules')
    [genes, rxnGeneMat] = getGenesFromGrRules(model.grRules);
    model.rxnGeneMat = sparse(numel(model.rxns),numel(model.genes));
    [~,geneOrder] = ismember(genes,model.genes);
    if any(geneOrder == 0)
        error(['The grRules includes the following gene(s), that are not in '...
            'the list of model genes: ', genes{~geneOrder}])
    end
    model.rxnGeneMat(:,geneOrder) = rxnGeneMat;
end

% Finalize GECKO model
if isGECKO
    % Fill in empty fields and empty entries
    for i={'kcat','source','notes','eccodes'} % Even keep empty
        model.ec = emptyOrFill(model.ec,i{1},{''},'rxns',true);
    end
    for i={'enzymes','mw','sequence'}
        model.ec = emptyOrFill(model.ec,i{1},{''},'genes',true);
    end
    model.ec = emptyOrFill(model.ec,'concs',{'NaN'},'genes',true);
    model.ec = emptyOrFill(model.ec,'kcat',{'0'},'genes',true);
    % Change string to double
    for i={'kcat','mw','concs'}
        if isfield(model.ec,i{1})
            model.ec.(i{1}) = str2double(model.ec.(i{1}));
        end
    end
    % Fill rxnEnzMat
    if ~isempty(enzStoich)
        rxnIdx              = cellfun('isempty', enzStoich(:,1));
        enzStoich(rxnIdx,:) = '';
        if ~isempty(enzStoich)
            rxnIdx              = cell2mat(enzStoich(:,1));
            [~,enzIdx]          = ismember(enzStoich(:,2),model.ec.enzymes);
            coeffs              = cell2mat(enzStoich(:,3));
            model.ec.rxnEnzMat  = zeros(numel(model.ec.rxns), numel(model.ec.genes));
            linearIndices       = sub2ind([numel(model.ec.rxns), numel(model.ec.genes)], rxnIdx, enzIdx);
            model.ec.rxnEnzMat(linearIndices) = coeffs;
        end
    end
    %Parse ec-codes (legacy multi-line eccodes path).
    if ~isempty(ecGecko)
        locs = cell2mat(ecGecko(:,1));
        for i=unique(locs)'
            ecGeckoCat=strjoin(ecGecko(locs==i,2),';');
            model.ec.eccodes{i,1}=ecGeckoCat;
        end
        emptyEc=cellfun('isempty',model.ec.eccodes);
        model.ec.eccodes(emptyEc)={''};
    end
end

if verbose
    fprintf(' Done!\n');
end
end

% --- Shared helpers ----------------------------------------------------

function model = emptyOrFill(model,field,emptyEntry,type,keepEmpty)
if nargin<5
    keepEmpty=false;
end
if isnumeric(emptyEntry)
    %A numeric field counts as "effectively empty" when it has no
    %elements at all, or when every element is NaN (matches the legacy
    %parser, which never allocated the field in the first place if no
    %row supplied a value).
    v = model.(field);
    emptyCells = isempty(v) || (isnumeric(v) && all(isnan(v(:))));
else
    emptyCells=cellfun('isempty',model.(field));
end
if all(emptyCells) && ~keepEmpty
    model = rmfield(model, field);
elseif numel(model.(field))<numel(model.(type))
    model.(field)(end+1:numel(model.(type)),1)=emptyEntry;
end
end
