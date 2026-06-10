function pairs = findDuplicateRxns(model, ignoreDirection)
% findDuplicateRxns  Find reactions that share identical stoichiometry.
%
% Counterpart of raven_python.manipulation.find_duplicate_reactions, and
% the upstream version of yeast-GEM's findDuplicatedRxns.
%
% Only stoichiometry is compared — bounds, GPRs, and annotations are
% ignored. The default treats A→B and B→A as duplicates (typical curation
% use case: "find reactions that could be merged").
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% ignoreDirection : logical, optional
%     treat A→B and B→A as duplicates (default true).
%
% Returns
% -------
% pairs : double
%     Nx2 numeric array of reaction-index pairs (i, j) where reactions i
%     and j share the same (possibly negated) stoichiometry, with i < j.
%     Empty if the model has no duplicates.
%
% Examples
% --------
%     pairs = findDuplicateRxns(model);
%     pairs = findDuplicateRxns(model, false);

if nargin < 2
    ignoreDirection = true;
end

pairs = zeros(0, 2);
n = numel(model.rxns);
for i = 1:n-1
    for j = i+1:n
        if isequal(model.S(:, i), model.S(:, j)) || ...
                (ignoreDirection && isequal(model.S(:, i), -model.S(:, j)))
            pairs(end+1, :) = [i, j]; %#ok<AGROW>
        end
    end
end
end
