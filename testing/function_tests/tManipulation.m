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

        function addRxnsCleanEquationDoesNotWarn(testCase)
            % The "metabolite on both sides" warning must only fire for the
            % equations it names, not on every call.
            rxnsToAdd.rxns      = {'newRxn'};
            rxnsToAdd.equations = {[testCase.model.mets{1} ' => ' testCase.model.mets{2}]};
            lastwarn('');
            evalc('addRxns(testCase.model, rxnsToAdd, 1);');
            testCase.verifyEmpty(strfind(lastwarn, 'both as substrate and product')); %#ok<STRIFY>
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
end
