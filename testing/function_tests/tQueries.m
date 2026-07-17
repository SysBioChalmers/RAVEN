classdef tQueries < RavenTestCase
% tQueries  Tests for the read-only query/accessor functions in queries/.

    methods (Test)

        function buildEquationReversible(testCase)
            eqn = buildEquation({'a';'b'}, [-1;1], true);
            testCase.verifyClass(eqn, 'char');
            testCase.verifySubstring(eqn, 'a');
            testCase.verifySubstring(eqn, 'b');
            testCase.verifySubstring(eqn, '<=>');
        end

        function buildEquationIrreversible(testCase)
            eqn = buildEquation({'a';'b'}, [-1;1], false);
            testCase.verifySubstring(eqn, '=>');
            testCase.verifyEmpty(strfind(eqn, '<=>')); %#ok<STRIFY>
        end

        function checkModelStructValidModel(testCase)
            % A valid model must not throw when errors are requested.
            % Advisory warnings (unused elements, bounds notes, etc.) are
            % expected MATLAB warnings in the new system; the check is
            % only that no error is raised.
            checkModelStruct(testCase.model, true);
        end

        function constructEquationsAllRxns(testCase)
            eqns = constructEquations(testCase.model);
            testCase.verifyClass(eqns, 'cell');
            testCase.verifyNumElements(eqns, numel(testCase.model.rxns));
            testCase.verifyTrue(all(~cellfun(@isempty, eqns)));
        end

        function constructEquationsSubset(testCase)
            eqns = constructEquations(testCase.model, testCase.model.rxns(1:3));
            testCase.verifyNumElements(eqns, 3);
        end

        function constructSSimple(testCase)
            [S, mets] = constructS({'a + b => c'});
            testCase.verifyEqual(sort(mets(:)), {'a';'b';'c'});
            testCase.verifySize(S, [3 1]);
            % reactants negative, product positive (S is sparse)
            testCase.verifyEqual(full(S), [-1;-1;1]);
        end

        function constructSLeadingNumberIsPartOfName(testCase)
            % A metabolite whose name starts with a number must not have that
            % number read as a stoichiometric coefficient when the metabolite
            % list says the whole entry is a metabolite.
            mets = {'2 oxoglutarate';'succinate'};
            [S, outMets] = constructS({'2 oxoglutarate => succinate'}, 'mets', mets);
            testCase.verifyEqual(outMets, mets);
            testCase.verifyEqual(full(S), [-1;1]);
        end

        function constructSLeadingNumberStillReadsCoefficient(testCase)
            % With no such metabolite, the leading number is a coefficient.
            mets = {'oxoglutarate';'succinate'};
            [S, outMets] = constructS({'2 oxoglutarate => succinate'}, 'mets', mets);
            testCase.verifyEqual(outMets, mets);
            testCase.verifyEqual(full(S), [-2;1]);
        end

        function getAllRxnsFromGenesType(testCase)
            % Use a reaction that has a gene association.
            withGpr = testCase.model.rxns(find(~cellfun(@isempty, ...
                testCase.model.grRules), 1));
            allRxns = getAllRxnsFromGenes(testCase.model, withGpr);
            testCase.verifyClass(allRxns, 'cell');
            testCase.verifyTrue(all(ismember(withGpr, allRxns)));
        end

        function getElementalBalanceFields(testCase)
            bs = getElementalBalance(testCase.model);
            testCase.verifyTrue(isfield(bs, 'balanceStatus'));
            testCase.verifyTrue(isfield(bs, 'leftComp'));
            testCase.verifyTrue(isfield(bs, 'rightComp'));
            testCase.verifySize(bs.balanceStatus, [numel(testCase.model.rxns) 1]);
        end

        function getElementalBalanceEmptyRxnIsUnbalanced(testCase)
            % A reaction with no metabolites (all-zero S column) must not
            % be falsely reported as balanced (B2 fix).
            m = testCase.model;
            m.rxns{end+1}    = 'emptyRxn';
            m.S(:, end+1)    = 0;
            m.lb(end+1)      = 0;
            m.ub(end+1)      = 1000;
            m.rev(end+1)     = 0;
            m.c(end+1)       = 0;
            m.grRules{end+1} = '';
            if isfield(m, 'rxnNames'),   m.rxnNames{end+1}  = 'emptyRxn'; end
            if isfield(m, 'rxnGeneMat'), m.rxnGeneMat(end+1,:) = 0;       end
            bs = getElementalBalance(m);
            testCase.verifyLessThanOrEqual(bs.balanceStatus(end), -1);
        end

        function getExchangeRxnsConsistent(testCase)
            [exch, idx] = getExchangeRxns(testCase.model);
            testCase.verifyClass(exch, 'cell');
            testCase.verifyNotEmpty(exch);
            testCase.verifyEqual(testCase.model.rxns(idx), exch);
        end

        function getExchangeRxnsThirdOutputIsMetIndex(testCase)
            [~, idx, mets] = getExchangeRxns(testCase.model);
            testCase.verifyNumElements(mets, numel(idx));
            % Each returned index must point to a metabolite with a
            % non-zero entry in the exchange reaction's S column.
            for i = 1:numel(idx)
                testCase.verifyNotEqual(mets(i), 0);
                testCase.verifyNotEqual(testCase.model.S(mets(i), idx(i)), 0);
            end
        end

        function getGenesFromGrRulesMatchesModel(testCase)
            [genes, rxnGeneMat] = getGenesFromGrRules(testCase.model.grRules);
            testCase.verifyTrue(iscellstr(genes)); %#ok<ISCLSTR>
            testCase.verifySize(rxnGeneMat, [numel(testCase.model.rxns) numel(genes)]);
            testCase.verifyEqual(sort(genes(:)), sort(testCase.model.genes(:)));
        end

        function getIndexesByName(testCase)
            idx = getIndexes(testCase.model, testCase.model.rxns(1:3), 'rxns');
            testCase.verifyEqual(idx, [1;2;3]);
        end

        function getIndexesLogical(testCase)
            lg = getIndexes(testCase.model, testCase.model.rxns(1), 'rxns', true);
            testCase.verifyClass(lg, 'logical');
            testCase.verifyEqual(find(lg), 1);
        end

        function getIndexesLogicalMaskAllTrueReturnsIndices(testCase)
            % A logical all-true mask of length n must return numeric [1..n],
            % not be passed through as a logical array (G5 fix: islogical check
            % instead of all() which cannot distinguish the two cases).
            nR = numel(testCase.model.rxns);
            idx = getIndexes(testCase.model, true(nR, 1), 'rxns');
            testCase.verifyClass(idx, 'double');
            testCase.verifyEqual(idx, (1:nR)');
        end

        function getRxnsInCompType(testCase)
            I = getRxnsInComp(testCase.model, 'c');
            testCase.verifyNotEmpty(I);
        end

        function getTransportRxnsType(testCase)
            tr = getTransportRxns(testCase.model);
            testCase.verifyClass(tr, 'logical');
            testCase.verifyNumElements(tr, numel(testCase.model.rxns));
        end

        function parseFormulasMW(testCase)
            [elements, useMat, exitFlag, MW] = parseFormulas(testCase.model.metFormulas); %#ok<ASGLU>
            testCase.verifyNumElements(MW, numel(testCase.model.mets));
            testCase.verifyTrue(isfield(elements, 'abbrevs'));
        end

        function parseRxnEquMetNames(testCase)
            mets = parseRxnEqu({'a + b => c'});
            testCase.verifyTrue(all(ismember({'a','b','c'}, mets)));
        end

        function modelSummaryNoFluxRuns(testCase)
            out = evalc('modelSummary(testCase.model)');
            testCase.verifyClass(out, 'char');
            testCase.verifySubstring(out, 'Reactions');
        end

        function modelSummaryWithFluxRuns(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            out = evalc('modelSummary(testCase.model, ''fluxes'', sol.x)');
            testCase.verifyClass(out, 'char');
            testCase.verifySubstring(out, 'Objective value');
        end

        function printFluxesRuns(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            out = evalc('printFluxes(testCase.model, sol.x, true)');
            testCase.verifyClass(out, 'char');
        end

        function printModelRuns(testCase)
            out = evalc('printModel(testCase.model, testCase.model.rxns(1))');
            testCase.verifyNotEmpty(out);
        end

        function printModelStatsRuns(testCase)
            out = evalc('printModelStats(testCase.model)');
            testCase.verifyClass(out, 'char');
        end

    end
end
