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
            % A valid model must not raise when errors are requested.
            testCase.verifyWarningFree(@() checkModelStruct(testCase.model, true));
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

        function getExchangeRxnsConsistent(testCase)
            [exch, idx] = getExchangeRxns(testCase.model);
            testCase.verifyClass(exch, 'cell');
            testCase.verifyNotEmpty(exch);
            testCase.verifyEqual(testCase.model.rxns(idx), exch);
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

        function getMetsInCompCount(testCase)
            I = getMetsInComp(testCase.model, 'c');
            testCase.verifyEqual(nnz(I), sum(testCase.model.metComps == 1));
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
