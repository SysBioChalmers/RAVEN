classdef tManipulation < RavenTestCase
% tManipulation  Tests for the structural model-editing functions in manipulation/.
%
%   Assertions favour structural invariants (counts, membership, round-trips,
%   idempotence) over full golden-model comparisons, so they stay robust while
%   still exercising each function's main path. evalc is used to suppress the
%   informational printing some functions do.

    methods (Test)

        function addExchangeRxnsAddsOne(testCase)
            [m2, added] = addExchangeRxns(testCase.model, 'both', '6pgl_c');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns) + 1);
            testCase.verifyNumElements(added, 1);
        end

        function addGenesRavenAddsGenes(testCase)
            g.genes = {'testgene1','testgene2'};
            g.geneShortNames = {'s1','s2'};
            m2 = addGenesRaven(testCase.model, g);
            testCase.verifyEqual(numel(m2.genes), numel(testCase.model.genes) + 2);
            testCase.verifyTrue(all(ismember(g.genes, m2.genes)));
        end

        function addMetsAddsMets(testCase)
            mta.metNames = {'newMetA','newMetB'};
            mta.compartments = {'c','e'};
            evalc('m2 = addMets(testCase.model, mta);');
            testCase.verifyEqual(numel(m2.mets), numel(testCase.model.mets) + 2);
        end

        function addRxnsAddsRxn(testCase)
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            evalc('m2 = addRxns(testCase.model, r, 2, ''c'', true);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns) + 1);
            testCase.verifyTrue(ismember('newRxn1', m2.rxns));
        end

        function addRxnsStringEqnTypeIdAlias(testCase)
            % 'id' is the string alias for eqnType=1 (match by model.mets).
            r.rxns = {'idAliasRxn'};
            r.equations = {[testCase.model.mets{1} ' => ' testCase.model.mets{2}]};
            m1 = addRxns(testCase.model, r, 1);
            m2 = addRxns(testCase.model, r, 'id');
            testCase.verifyEqual(m1.S, m2.S);
        end

        function addRxnsStringEqnTypeNameAlias(testCase)
            % 'name' is the string alias for eqnType=2 (match by model.metNames).
            r.rxns = {'nameAliasRxn'};
            r.equations = {'2-Oxoglutarate => NAMETEST'};
            evalc('m1 = addRxns(testCase.model, r, 2, ''c'', true);');
            evalc('m2 = addRxns(testCase.model, r, ''name'', ''c'', true);');
            testCase.verifyEqual(m1.S, m2.S);
        end

        function addRxnsStringEqnTypeNameCompAlias(testCase)
            % 'name[comp]' is the string alias for eqnType=3.
            r.rxns = {'namecompAliasRxn'};
            r.equations = {'2-Oxoglutarate[c] => COMPTEST[c]'};
            evalc('m1 = addRxns(testCase.model, r, 3, [], true);');
            evalc('m2 = addRxns(testCase.model, r, ''name[comp]'', [], true);');
            testCase.verifyEqual(m1.S, m2.S);
        end

        function addRxnsInvalidStringEqnTypeErrors(testCase)
            % An unrecognised string for eqnType must throw.
            r.rxns = {'errRxn'};
            r.equations = {[testCase.model.mets{1} ' => ' testCase.model.mets{2}]};
            testCase.verifyError(@() addRxns(testCase.model, r, 'invalid'), ?MException);
        end

        function addRxnsGenesMetsCopiesFromSource(testCase)
            sbmlFile = fullfile(testCase.ravenRoot,'tutorial','empty.xml');
            testCase.assumeDependency(exist(sbmlFile,'file')==2, 'tutorial/empty.xml');
            evalc('src = importModel(sbmlFile, [], true);');
            evalc('m2 = addRxnsGenesMets(testCase.model, src, ''r1'', true);');
            testCase.verifyTrue(ismember('r1', m2.rxns));
        end

        function addTransportAddsRxn(testCase)
            evalc(['m2 = addTransport(testCase.model, ''c'', ''e'', ' ...
                '{''6-phospho-D-glucono-1,5-lactone''}, false, false, ''tr_'');']);
            testCase.verifyGreaterThan(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function changeGrRulesUpdatesRule(testCase)
            m2 = changeGrRules(testCase.model, 'ACKr', 'b2296 and b1849', true);
            idx = strcmp(m2.rxns, 'ACKr');
            testCase.verifyEqual(m2.grRules{idx}, 'b2296 and b1849');
        end

        function changeRxnsUpdatesEquation(testCase)
            evalc('m2 = changeRxns(testCase.model, ''ACKr'', ''2-Oxoglutarate <=> TEST'', 2, ''c'', true);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function closeModelReturnsStruct(testCase)
            m2 = closeModel(testCase.model);
            testCase.verifyClass(m2, 'struct');
            testCase.verifyGreaterThanOrEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function contractModelNoMoreRxns(testCase)
            evalc('m2 = contractModel(testCase.model);');
            testCase.verifyLessThanOrEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function convertToIrrevAllIrreversible(testCase)
            m2 = convertToIrrev(testCase.model);
            testCase.verifyTrue(all(m2.rev == 0));
            testCase.verifyGreaterThanOrEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function convertToIrrevSplitsBoundsAndStoichiometry(testCase)
            % A reversible reaction with bounds (-500,1000) keeps (0,1000) and
            % its stoichiometry, while the _REV copy gets (0,500) and the
            % negated stoichiometry.
            m = tManipulation.twoMetModel();
            m.rev = 1; m.lb = -500; m.ub = 1000;
            m2 = convertToIrrev(m);

            fwd = strcmp(m2.rxns, 'R1');
            rev = strcmp(m2.rxns, 'R1_REV');
            testCase.verifyTrue(any(rev));
            testCase.verifyEqual(full(m2.lb(fwd)), 0);
            testCase.verifyEqual(full(m2.ub(fwd)), 1000);
            testCase.verifyEqual(full(m2.lb(rev)), 0);
            testCase.verifyEqual(full(m2.ub(rev)), 500);
            testCase.verifyEqual(full(m2.S(:,fwd)), [-1;1]);
            testCase.verifyEqual(full(m2.S(:,rev)), [1;-1]);
        end

        function convertToIrrevReverseCopyInheritsGrRule(testCase)
            % The _REV copy is catalysed by the same genes as the forward one.
            m = tManipulation.twoMetModel();
            m.rev = 1; m.lb = -500; m.ub = 1000;
            m.genes = {'g1'}; m.grRules = {'g1'}; m.rxnGeneMat = sparse(1,1,1);
            m2 = convertToIrrev(m);
            testCase.verifyEqual(m2.grRules{strcmp(m2.rxns,'R1_REV')}, 'g1');
        end

        function findDuplicateRxnsIgnoreDirection(testCase)
            % a -> b and b -> a are the same reaction run backwards, so they
            % group by default and stay separate when direction matters.
            m = tManipulation.twoMetModel();
            m.rxns     = {'R1';'R2'};
            m.rxnNames = {'R1';'R2'};
            m.S   = sparse([-1 1; 1 -1]);
            m.lb  = [0;0]; m.ub = [1000;1000]; m.rev = [0;0]; m.c = [0;0];
            m.grRules = {'';''}; m.rxnGeneMat = sparse(2,0);

            pairs = findDuplicateRxns(m);
            testCase.verifyEqual(pairs, [1 2]);

            pairs = findDuplicateRxns(m, 'ignoreDirection', false);
            testCase.verifyEmpty(pairs);
        end

        function changeGrRulesAppendsToExistingRule(testCase)
            % replace=false must OR the new rule onto the existing one rather
            % than overwrite it, and add the new gene to the model.
            m2 = changeGrRules(testCase.model, 'ACKr', 'b9999', false);
            rule = m2.grRules{strcmp(m2.rxns, 'ACKr')};
            testCase.verifySubstring(rule, 'b9999');
            testCase.verifySubstring(rule, ' or ');
            % the gene the reaction already had must still be in the rule
            oldRule = testCase.model.grRules{strcmp(testCase.model.rxns, 'ACKr')};
            oldGene = regexp(oldRule, 'b\d+', 'match', 'once');
            testCase.verifySubstring(rule, oldGene);
            testCase.verifyTrue(ismember('b9999', m2.genes));
        end

        function copyToCompsDeleteOriginalIsAMove(testCase)
            % deleteOriginal turns the copy into a move: the reaction count is
            % unchanged and the original compartment's copy is gone.
            evalc('copied = copyToComps(testCase.model, {''p''}, ''rxns'', ''ACKr'');');
            evalc(['moved = copyToComps(testCase.model, {''p''}, ''rxns'', ''ACKr'', ' ...
                '''deleteOriginal'', true);']);
            testCase.verifyEqual(numel(copied.rxns), numel(testCase.model.rxns) + 1);
            testCase.verifyEqual(numel(moved.rxns), numel(testCase.model.rxns));
        end

        function mergeModelsMetParamDecidesUnification(testCase)
            % The same metabolite under two ids unifies when matching on names
            % and stays distinct when matching on ids.
            a = tManipulation.namedMetModel('glc_c', 'A');
            b = tManipulation.namedMetModel('glucose_c', 'B');
            evalc('byName = mergeModels({a; b});');
            evalc('byId   = mergeModels({a; b}, ''metParam'', ''mets'');');
            testCase.verifyEqual(nnz(strcmp(byName.metNames, 'Glucose')), 1);
            testCase.verifyEqual(nnz(strcmp(byId.metNames, 'Glucose')), 2);
        end

        function copyToCompsAddsCompartment(testCase)
            evalc('m2 = copyToComps(testCase.model, {''p''}, ''ACKr'');');
            testCase.verifyTrue(ismember('p', m2.comps));
            testCase.verifyGreaterThan(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function deleteUnusedGenesNoMoreGenes(testCase)
            evalc('m2 = deleteUnusedGenes(testCase.model);');
            testCase.verifyLessThanOrEqual(numel(m2.genes), numel(testCase.model.genes));
        end

        function expandModelNoFewerRxns(testCase)
            evalc('m2 = expandModel(testCase.model);');
            testCase.verifyGreaterThanOrEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function findDuplicateRxnsReturnsResult(testCase)
            pairs = findDuplicateRxns(testCase.model);
            testCase.verifyTrue(isnumeric(pairs) || islogical(pairs) || iscell(pairs));
        end

        function generateNewIdsUniqueAndPrefixed(testCase)
            ids = generateNewIds(testCase.model, 'rxns', 'r_', 5);
            testCase.verifyNumElements(ids, 5);
            testCase.verifyEqual(numel(unique(ids)), 5);
            testCase.verifyTrue(all(startsWith(ids, 'r_')));
        end

        function mergeCompartmentsSingleComp(testCase)
            evalc('m2 = mergeCompartments(testCase.model);');
            testCase.verifyNumElements(m2.comps, 1);
        end

        function mergeModelsReturnsStruct(testCase)
            modelB = removeReactions(testCase.model, testCase.model.rxns(1:5), true, true, true);
            evalc('merged = mergeModels({testCase.model; modelB});');
            testCase.verifyClass(merged, 'struct');
            testCase.verifyNotEmpty(merged.rxns);
        end

        function permuteModelReversesOrder(testCase)
            n = numel(testCase.model.rxns);
            m2 = permuteModel(testCase.model, n:-1:1, 'rxns');
            testCase.verifyEqual(m2.rxns, flipud(testCase.model.rxns));
        end

        function permuteModelRemapsRxnComps(testCase)
            % After swapping compartments the rxnComps and geneComps indices
            % in the OUTPUT must follow the permutation; the input is unchanged.
            m = struct();
            m.rxns = {'r1'; 'r2'};
            m.mets = {'a_c'; 'b_e'};
            m.S    = sparse([1 -1; -1 1]);
            m.lb   = [0; 0]; m.ub = [1000; 1000]; m.rev = [0; 0];
            m.c    = [1; 0]; m.b  = [0; 0];
            m.comps     = {'c'; 'e'};
            m.compNames = {'Cytoplasm'; 'Extracellular'};
            m.metComps  = [1; 2];
            m.genes     = {'g1'; 'g2'};
            m.rxnGeneMat = sparse(eye(2));
            m.grRules   = {'g1'; 'g2'};
            m.rxnComps  = [1; 2];   % r1→comps{1}, r2→comps{2}
            m.geneComps = [1; 2];
            m2 = permuteModel(m, [2, 1], 'comps');
            % Swap: old comps{1}→new position 2, old comps{2}→new position 1
            testCase.verifyEqual(m2.rxnComps,  [2; 1]);
            testCase.verifyEqual(m2.geneComps, [2; 1]);
            % Input must be unmodified
            testCase.verifyEqual(m.rxnComps,  [1; 2]);
            testCase.verifyEqual(m.geneComps, [1; 2]);
        end

        function removeBadRxnsReturnsStruct(testCase)
            evalc('m2 = removeBadRxns(testCase.model);');
            testCase.verifyClass(m2, 'struct');
        end

        function removeGenesRemovesGene(testCase)
            m2 = removeGenes(testCase.model, 'b1817', true, true, false);
            testCase.verifyFalse(ismember('b1817', m2.genes));
            testCase.verifyLessThan(numel(m2.genes), numel(testCase.model.genes));
        end

        function removeMetsRemovesMet(testCase)
            m2 = removeMets(testCase.model, 'Acetate', true, true, true, true);
            testCase.verifyLessThan(numel(m2.mets), numel(testCase.model.mets));
        end

        function removeReactionsRemovesRxn(testCase)
            m2 = removeReactions(testCase.model, testCase.model.rxns(1), true, true, true);
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns) - 1);
            testCase.verifyFalse(ismember(testCase.model.rxns{1}, m2.rxns));
        end

        function replaceMetsMergesMet(testCase)
            evalc('m2 = replaceMets(testCase.model, ''Acetaldehyde'', ''Acetate'');');
            testCase.verifyLessThan(numel(m2.mets), numel(testCase.model.mets));
        end

        function setExchangeBoundsRuns(testCase)
            evalc('m2 = setExchangeBounds(testCase.model, {''ac_e'';''akg_e''}, -500, 500);');
            testCase.verifyClass(m2, 'struct');
        end

        function setParamObjective(testCase)
            m2 = setParam(testCase.model, 'obj', testCase.model.rxns(1), 1);
            idx = strcmp(m2.rxns, testCase.model.rxns{1});
            testCase.verifyEqual(m2.c(idx), 1);
        end

        function setParamUpperBound(testCase)
            m2 = setParam(testCase.model, 'ub', testCase.model.rxns(1), 42);
            idx = strcmp(m2.rxns, testCase.model.rxns{1});
            testCase.verifyEqual(m2.ub(idx), 42);
        end

        function simplifyModelReturnsStruct(testCase)
            evalc('m2 = simplifyModel(testCase.model);');
            testCase.verifyClass(m2, 'struct');
            testCase.verifyLessThanOrEqual(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function sortIdentifiersSortsRxns(testCase)
            m2 = sortIdentifiers(testCase.model);
            testCase.verifyEqual(sort(m2.rxns), m2.rxns);
        end

        function sortModelPreservesCounts(testCase)
            m2 = sortModel(testCase.model);
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns));
            testCase.verifyEqual(numel(m2.mets), numel(testCase.model.mets));
        end

        function standardizeGrRulesReturnsRules(testCase)
            evalc('grRules = standardizeGrRules(testCase.model);');
            testCase.verifyNumElements(grRules, numel(testCase.model.rxns));
        end

    end

    methods (Static, Access = private)

        function m = twoMetModel()
            % Single reaction a -> b in one compartment.
            m = struct();
            m.id        = 'toy';
            m.rxns      = {'R1'};
            m.rxnNames  = {'R1'};
            m.mets      = {'a';'b'};
            m.metNames  = {'a';'b'};
            m.metComps  = [1;1];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-1;1]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0; m.b = zeros(2,1);
            m.genes = {}; m.grRules = {''}; m.rxnGeneMat = sparse(1,0);
            m.metFormulas = {'C';'C'};
        end

        function m = namedMetModel(glucoseId, modelId)
            % Glucose[c] under a caller-chosen id, consumed by one reaction.
            m = struct();
            m.id        = modelId;
            m.rxns      = {['R_' modelId]};
            m.rxnNames  = m.rxns;
            m.mets      = {glucoseId; ['product_' modelId]};
            m.metNames  = {'Glucose'; ['Product' modelId]};
            m.metComps  = [1;1];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-1;1]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0; m.b = zeros(2,1);
            m.genes = {}; m.grRules = {''}; m.rxnGeneMat = sparse(1,0);
            m.metFormulas = {'C6H12O6';'C'};
        end

    end
end
