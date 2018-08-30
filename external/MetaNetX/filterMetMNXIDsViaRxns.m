function [fmodel,removed] = filtermetMetaNetXIDsViaRxns(model,mnx,ignoreComp,keepAtLeastOne)
%filtermetMetaNetXIDsViaRxns  Remove met MetaNetXID associations based on their rxns.
%
% filtermetMetaNetXIDsViaRxns determines which met MetaNetXIDs (if any) should be
% removed for mets associated with multiple MetaNetXIDs. This is accomplished by
% finding all model rxns that involve the met, obtaining their associated
% MNX rxn IDs, and retrieving those rxns from the MNX database. Any MetaNetXIDs
% associated with the met that are not included in that set of rxns from
% the MNX database will be removed.
%
% NOTE: metabolites that occur only in reactions that have no associated
%       rxnMetaNetXIDs will be skipped, as there is insufficient evidence to
%       properly filter their associated metMetaNetXIDs.
%
%
% USAGE:
%
%   [fmodel,removed] = filtermetMetaNetXIDsViaRxns(model,mnx,ignoreComp,keepAtLeastOne);
%
% INPUT:
%
%   model    A model structure containing reaction and metabolite MNX ID
%            association fields ("rxnMetaNetXID" and "metMetaNetXID", respectively).
%
%   mnx      (Optional) An MNX database structure, containing reaction-
%            related information retrieved from the MNX database, generated
%            using the following command: mnx = buildMNXref('rxns');
%            By default, the function will automatically run the above
%            command to regenerate the MNX database structure (slower).
%
%   ignoreComp   (Optional, Default FALSE) If TRUE, metabolite compartments
%                will be ignored. In this case, identical mets of different
%                compartments will be lumped together, and when searching
%                for rxns involving the metabolite, the compartment will be
%                ignored.
%
%   keepAtLeastOne  (Optional, Default FALSE) In some cases, none of the
%                   MetaNetXIDs associated with a metabolite are found in any of
%                   the reactions, and will result in a met with zero MetaNetXID
%                   associations in the filtered model.
%                   If keepAtLeastOne = TRUE, then in cases such as this,
%                   the MetaNetXID associations will not be removed in the
%                   filtered model, but will be indicated in the "removed"
%                   structure, with the text "ALL SHOULD BE REMOVED!" next
%                   to the metabolite ID.
%
% OUTPUT:
%
%   fmodel   A filtered model, which has removed metabolite MetaNetXIDs that
%            were not found in any of the MNX reactions associated with the
%            model and that metabolite.
%
%   removed  A structure containing more detailed information on the met
%            MetaNetXID associations that were removed.
%
%
% Jonathan Robinson 2018-05-28


% handle input arguments
if nargin < 2 || isempty(mnx)
    mnx = buildMNXref('rxns');
end
if nargin < 3
    ignoreComp = false;
end
if nargin < 4
    keepAtLeastOne = false;
end

% initialize "removed" structure and "fmodel" outputs
fmodel = model;
removed.mets = {};
removed.metMetaNetXID = {};

% if metMetaNetXID and/or rxnMetaNetXID field contains multiple columns, consolidate
% into a single column with nested cell entries
if size(model.metMetaNetXID,2) > 1
    model.metMetaNetXID = nestCell(model.metMetaNetXID,true);
end
if size(model.rxnMetaNetXID,2) > 1
    model.rxnMetaNetXID = nestCell(model.rxnMetaNetXID,true);
end

if ( ignoreComp )
    
    % check if model has already merged met compartments, or if "mets"
    % contains non-unique elements
    if length(unique(model.mets)) ~= length(model.mets)
        error('Input model "mets" field must contain unique (non-repeated) elements.');
    elseif any(regexp(model.mets{1},'\d$'))
        % Note: this test is specific to models whose met IDs end in
        % numbers (e.g., "m00001"), with compartments appended to the end
        % of the ID (e.g., "m00001c" or "m00001[c]").
        fprintf('\nIt appears that the model compartments have already been merged.\n');
        fprintf('The "ignoreComp" flag will be ignored, since it will have no effect.\n');
        S = model.S;
        ignoreComp = false;
    else
        
        % strip compartment label from model met ID
        if endsWith(model.mets{1},']')
            % compartment name is formatted as "m00001[c]"
            model.mets = regexprep(model.mets,'\[\w\]$','');
        else
            % compartment name is formatted as "m00001c"
            model.mets = regexprep(model.mets,'.$','');
        end
        
        % check if ignoring compartments will actually do anything
        if length(unique(model.mets)) == length(model.mets)
            fprintf('\nIt appears that the model compartments have already been merged.\n');
            fprintf('The "ignoreComp" flag will be ignored, since it will have no effect.\n');
            S = model.S;
            ignoreComp = false;
        else
            
            fprintf('\nMetabolite compartments will be ignored. NOTE: This process will merge metMetaNetXIDs\n');
            fprintf('of mets that are identical except for their compartment. If any mets are associated\n');
            fprintf('with compartment-specific metMetaNetXIDs, this is not recommended.\n');
            
            % If ignoring compartments, convert the stoich matrix into a
            % binary met-rxn association matrix (i.e., set all nonzero
            % entries = 1), and combine all associations for metabolites
            % that are identical except for their compartment. Also combine
            % their metMetaNetXIDs.
            S = (model.S ~= 0);
            [~,uniq_ind,met_groups] = unique(model.mets);
            h = waitbar(0,'Merging model metabolites across compartments...');
            for i = 1:max(met_groups)
                ind = (met_groups == i);
                S(ind,:) = repmat(any(S(ind,:),1),sum(ind),1);
                model.metMetaNetXID(ind) = repmat({unique(horzcat(model.metMetaNetXID{ind}))},sum(ind),1);
                waitbar(i/max(met_groups),h);
            end
            close(h);
            
            % now merge mets of different compartments into single met
            model.mets_nocomp = model.mets;  % first save these for indexing later
            model.mets = model.mets(uniq_ind);
            S = S(uniq_ind,:);
            model.metMetaNetXID = model.metMetaNetXID(uniq_ind);
            
        end
        
    end
    
else
    S = model.S;
end

% % Identify all mets with two or more MetaNetXID associations. Mets with one or
% % zero MetaNetXIDs will be ignored.
% multi_ind = find(cellfun(@numel,model.metMetaNetXID) > 1);

% Iterate through mets with multiple MetaNetXIDs, and determine which MetaNetXIDs
% should be removed.
h = waitbar(0,'Processing metabolites...');
for i = 1:length(model.mets)
    
    % get indices of all rxns in which the met participates
    rxn_ind = S(i,:) ~= 0;
    
    % obtain a unique list of all the rxn MetaNetXIDs associated with those
    % rxns, as well as their indices in the MNX database structure
    MNX_rxnIDs = unique(vertcat(model.rxnMetaNetXID{rxn_ind}));
    MNX_rxn_ind = ismember(mnx.rxns,MNX_rxnIDs);
    
    if isempty(MNX_rxnIDs)
        % If none of the rxns are associated to any rxn MetaNetXIDs, just skip
        % this metabolite. Otherwise, all met MetaNetXIDs associated with the
        % metabolite will be removed, which is not so helpful.
        continue
    end
    
    % obtain unique set of all mets participating in the MNX rxns
    MNX_metIDs = unique([mnx.rxnMets{MNX_rxn_ind}]);
    
    % determine which of the MetaNetXIDs currently associated to the model met
    % are not found in any of the MNX rxns
    rem_mets = ~ismember(model.metMetaNetXID{i},MNX_metIDs);
    
    if ~any(rem_mets)
        % if no met MetaNetXID associations should be removed, continue to next met
        continue
    elseif all(rem_mets) && (keepAtLeastOne)
        % In this case, all the MetaNetXIDs associated with this met are to be
        % removed, but the user has indicated that they do not want this to
        % happen. Therefore,  NONE will be removed, but the information
        % will be reported in the "removed" structure, so these cases can
        % be manually evaluated by the user.
        removed.mets = [removed.mets; model.mets(i)];
        removed.metMetaNetXID = [removed.metMetaNetXID; {'ALL SHOULD BE REMOVED!'}];
    else
        % remove one or more MetaNetXID associations from the metabolite
        removed.mets = [removed.mets; model.mets(i)];
        removed.metMetaNetXID = [removed.metMetaNetXID; {model.metMetaNetXID{i}(rem_mets)}];
        model.metMetaNetXID{i}(rem_mets) = [];
    end
    
    waitbar(i/length(model.mets),h);
end
close(h);

if any(cellfun(@(x) ismember({'ALL SHOULD BE REMOVED!'},x),removed.metMetaNetXID))
    fprintf('\n*** NOTE: There was at least one case where ALL the MetaNetXIDs associated\n');
    fprintf('          with a metabolite were not found in any of the involved rxns,\n');
    fprintf('          and therefore NONE were removed. See "removed" structure for\n');
    fprintf('          further information.\n\n');
end


% update filtered model (fmodel) metMetaNetXID field
if ( ignoreComp )
    % if compartment was ignored (and thus mets were merged), distribute
    % the updated metMetaNetXID assignments to all compartment-versions of each
    % metabolite.
    for i = 1:length(model.mets)
        ind = ismember(model.mets_nocomp,model.mets(i));
        fmodel.metMetaNetXID(ind) = model.metMetaNetXID(i);
    end
else
    fmodel.metMetaNetXID = model.metMetaNetXID;
end
