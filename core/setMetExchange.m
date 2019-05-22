function [exchModel,unusedMets] = setMetExchange(model,mets,lb,ub,closeOthers,mediaOnly)
% setMetExchange
%   Define the exchange flux bounds for a set of metabolites.
%
% Input:
%   model         a model structure
%   mets          a cell array of metabolite names (case insensitive) or 
%                 metabolite IDs, or a vector of metabolite indices
%                 (opt, default all exchanged metabolites)
%   lb            lower bound of exchange flux. Can be either a vector of
%                 bounds corresponding to each of the provided metabolites,
%                 or a single value that will be applied to all.
%                 (opt, default to model.annotation.defaultLB if it exists,
%                 otherwise -1000)
%   ub            upper bound of exchange flux. Can be either a vector of
%                 bounds corresponding to each of the provided metabolites,
%                 or a single value that will be applied to all.
%                 (opt, default to model.annotation.defaultLB if it exists,
%                 otherwise 1000)
%   closeOthers   close exchange reactions for all other exchanged 
%                 metabolites not present in the provided list. This will
%                 prevent IMPORT of the metabolites, but their EXPORT will
%                 not be modified.
%                 (opt, default true)
%   mediaOnly     only consider exchange reactions involving exchange to or
%                 from the extracellular (media) compartment. Reactions
%                 such as "sink" reactions that exchange metabolites
%                 directly with an intracellular compartment will therefore
%                 be ignored even though "getExchangeRxns" identifies such
%                 such reactions as exchange reactions.
%                 Note: The function will attempt to identify the
%                 extracellular compartment by the "compNames" field, and
%                 also requires the "metComps" field to be present,
%                 otherwise the mediaOnly flag will be ignored.
%                 (opt, default false)
%
% Output:
%   exchModel     a model structure with updated exchange flux bounds for
%                 the provided set of metabolites
%   unusedMets    metabolites provided by the user that were not used
%                 because they are not involved in any exchange reactions
%                 in the model
%
% NOTE: Exchange reactions involving more than one metabolite will be
% ignored.
%
% Usage: exchModel = setMetExchange(model,mets,lb,ub,closeOthers,mediaOnly);
%
%
% Jonathan Robinson, 2019-05-22
%


% handle input arguments
if nargin < 2
    mets = [];
elseif ischar(mets)
    mets = {mets};  % in case only one metabolite is provided as a string
end 

if nargin < 3 || isempty(lb)
    if isfield(model,'annotation') && isfield(model.annotation,'defaultLB')
        lb = model.annotation.defaultLB;
    else
        lb = -1000;
    end
end

if nargin < 4 || isempty(ub)
    if isfield(model,'annotation') && isfield(model.annotation,'defaultUB')
        ub = model.annotation.defaultUB;
    else
        ub = 1000;
    end
end

if nargin < 5 || isempty(closeOthers)
    closeOthers = true;
end

if nargin < 6
    mediaOnly = false;
elseif mediaOnly
    if ~all(isfield(model,{'compNames','metComps'}))
        error('mediaOnly option requires the "compNames" and "metComps" model fields.');
    end
end


% for models with an "unconstrained" field, generate a version of the model
% where all the unconstrained metabolites have stoichiometric coeffs set to
% zero, but remain in the model to avoid index changes
model_temp = model;
if isfield(model,'unconstrained')
    model_temp.S(model_temp.unconstrained == 1,:) = 0;
end

% find exchange rxns, ignoring those involving more than one metabolite
[~,exchRxnInd] = getExchangeRxns(model);
exchRxnInd(sum(model_temp.S(:,exchRxnInd) ~= 0,1) > 1) = [];

% find all exchanged metabolites
exchMetInd = find(any(model_temp.S(:,exchRxnInd) ~= 0,2));

if ( mediaOnly )
    % ignore exchanged metabolites in non-extracellular compartments and
    % any exchange reactions involving these metabolites
    [~,extCompInd] = ismember('extracellular',lower(model.compNames));
    if extCompInd > 0
        ignoreMet = (model.metComps(exchMetInd) ~= extCompInd);
    else
        error('Could not find any compartments named "extracellular".');
    end
    ignoreRxn = any(model.S(exchMetInd(ignoreMet),exchRxnInd) ~= 0,1);
    exchMetInd(ignoreMet) = [];
    exchRxnInd(ignoreRxn) = [];
end

% Check that all exchange reactions are formulated in the same direction.
% If not, this means that negative flux indicates import for some exchange
% reactions, but indicates export for others. Therefore, the LB and UB
% would need to be specified differently depending on the exchange reaction
% direction, which is error-prone.
if all(all(model.S(exchMetInd,exchRxnInd) <= 0))
    importDir = 'backward';
elseif all(all(model.S(exchMetInd,exchRxnInd) >= 0))
    importDir = 'forward';
else
    fprintf('WARNING: Some exchange reactions differ in direction, and therefore have opposite meanings of LB and UB.');
    if closeOthers
        fprintf('         Therefore, the "closeOthers" option will be set to FALSE.');
    end
end

% prepare exchanged metabolites and bounds
if ~isempty(mets)
    
    % prepare lb and ub variables
    if numel(lb) == 1
        lb = lb*ones(size(exchMetInd));
    elseif numel(lb) ~= numel(mets)
        error('lb must be a single value or a vector of equal length as mets');
    end
    if numel(ub) == 1
        ub = ub*ones(size(exchMetInd));
    elseif numel(ub) ~= numel(mets)
        error('ub must be a single value or a vector of equal length as mets');
    end
    
    % map provided mets to exchanged metabolites
    if isnumeric(mets)
        % mets are provided as indices
        exchMetInd = intersect(mets,exchMetInd);
        [~,keepMet] = ismember(exchMetInd,mets);
    elseif sum(ismember(lower(mets),lower(model.metNames))) > sum(ismember(mets,model.mets))
        % assume that mets are provided as names
        metInd = find(ismember(lower(model.metNames),lower(mets)));
        exchMetInd = intersect(metInd,exchMetInd);
        [~,keepMet] = ismember(lower(model.metNames(exchMetInd)),lower(mets));
    else
        % assume that mets are provided as met IDs
        metInd = find(ismember(model.mets,mets));
        exchMetInd = intersect(metInd,exchMetInd);
        [~,keepMet] = ismember(model.mets(exchMetInd),mets);
    end
    
    % update bound vectors
    lb = lb(keepMet);
    ub = ub(keepMet);
    
    % get provided metabolites not used
    unusedMets = setdiff(mets,mets(keepMet));
    
else
    % if no mets were provided, use all exchanged mets
    if numel(lb) > 1 || numel(ub) > 1
        error('Only one upper and one lower bound may be provide if metabolites are not specified.');
    else
        lb = lb*ones(size(exchMetInd));
        ub = ub*ones(size(exchMetInd));
    end
    unusedMets = {};
end

% check that at least one exchanged metabolite matches
if isempty(exchMetInd)
    fprintf('None of the provided metabolites were found in any exchange reactions.\n');
    exchModel = model;
    return
end

% determine which metabolite is exchanged in each exchange reaction
[metInd,rxnInd] = find(model_temp.S(exchMetInd,exchRxnInd) ~= 0);

% check for any metabolites that are exchanged in more than one reaction
tbl = tabulate(metInd);
repeatedInds = tbl(:,2) > 1;
multiMetInd = exchMetInd(metInd(repeatedInds));
if ~isempty(multiMetInd)
    fprintf('WARNING: The following metabolites are involved in more than one exchange reaction:\n');
    fprintf('\t%s\n',model.metNames{multiMetInd(1:min(numel(multiMetInd),10))});
    if numel(multiMetInd) > 10
        fprintf('\t... and %u more.\n',numel(multiMetInd)-10);
    end
end

% set exchange reaction bounds
model = setParam(model,'lb',exchRxnInd(rxnInd),lb(metInd));
model = setParam(model,'ub',exchRxnInd(rxnInd),ub(metInd));

if closeOthers
    % constrain import of all other exchange reactions to zero
    constrainInd = setdiff(exchRxnInd,exchRxnInd(rxnInd));
    if strcmp(importDir,'backward')
        model = setParam(model,'lb',constrainInd,0);
    else
        model = setParam(model,'ub',constrainInd,0);
    end
end

% assign output model
exchModel = model;




