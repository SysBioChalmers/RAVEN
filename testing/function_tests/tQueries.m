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

        function checkModelStructNoFalsePositiveOnWordStartingName(testCase)
            % A metabolite name beginning with a non-numeric word must not be
            % flagged as "begins with a number": str2double of that word is
            % NaN, and any(NaN) is true in MATLAB.
            m = testCase.model;
            m.metNames{1} = 'alpha keto acid';
            issues = checkModelStruct(m);
            hit = arrayfun(@(x) contains(x.message,'begin with a number'), issues);
            testCase.verifyFalse(any(hit));
        end

        function checkModelStructFlagsDuplicateNonSboMiriam(testCase)
            % Two metabolites with different names sharing the same
            % non-SBO MIRIAM (e.g. the same KEGG id) must be flagged:
            % regexp's [] "no match" result must not be mistaken for a
            % match via bitwise negation (~[] is also [], so the check
            % never fired at all before the fix).
            m = testCase.model;
            idx = find(~strcmp(m.metNames, m.metNames{1}), 1);
            m.metMiriams = cell(numel(m.mets),1);
            m.metMiriams{1}.name = {'kegg.compound'}; m.metMiriams{1}.value = {'C00031'};
            m.metMiriams{idx}.name = {'kegg.compound'}; m.metMiriams{idx}.value = {'C00031'};
            issues = checkModelStruct(m);
            hit = arrayfun(@(x) contains(x.message,'more than one unique metabolite name'), issues);
            testCase.verifyTrue(any(hit));
        end

        function checkModelStructThrowErrorsFalseSkipsAfterMissingField(testCase)
            % throwErrors=false must warn about a missing required field
            % without then crashing trying to read that very field in a
            % later check.
            m = rmfield(testCase.model, 'id');
            evalc('checkModelStruct(m, ''throwErrors'', false);');
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

        function constructSSumsRepeatedMetabolite(testCase)
            % A metabolite written more than once on the same side is one
            % coefficient, not the last occurrence.
            [S, mets] = constructS({'2 a + a => b'}, 'mets', {'a';'b'});
            testCase.verifyEqual(mets, {'a';'b'});
            testCase.verifyEqual(full(S), [-3;1]);
        end

        function constructSCancelsMetaboliteOnBothSides(testCase)
            % A metabolite on both sides cancels, which is exactly what
            % badRxns reports: the reaction is empty in the S matrix.
            [S, ~, badRxns] = constructS({'atp + adp <=> adp + atp'}, ...
                'mets', {'atp';'adp'});
            testCase.verifyEqual(full(S), [0;0]);
            testCase.verifyTrue(badRxns(1));
        end

        function constructSMissingMetKeepsPercentInErrorNoRxns(testCase)
            % A missing metabolite name containing "%" must survive intact
            % in the error message, not be truncated by sprintf misreading
            % it as a format directive.
            testCase.verifyError(@() constructS({'a + met_with_%_sign => b'}, ...
                'mets', {'a';'b'}), 'RAVEN:badInput');
            try
                constructS({'a + met_with_%_sign => b'}, 'mets', {'a';'b'});
            catch e
                testCase.verifySubstring(e.message, 'met_with_%_sign');
            end
        end

        function constructSMissingMetKeepsPercentInErrorWithRxns(testCase)
            % Same, but via the reaction-annotated branch (rxns supplied),
            % which also splices in the reaction id.
            try
                constructS({'a + met_with_%_sign => b'}, 'mets', {'a';'b'}, ...
                    'rxns', {'rxn_%_id'});
                testCase.verifyFail('Expected an error to be thrown.');
            catch e
                testCase.verifySubstring(e.message, 'met_with_%_sign');
                testCase.verifySubstring(e.message, 'rxn_%_id');
            end
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
            % be falsely reported as balanced.
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

        function getElementalBalanceReportsChargeBalance(testCase)
            % A charge-imbalanced reaction is reported as such, while the
            % elemental verdict is left alone.
            m = testCase.model;
            m.mets      = {'a';'b'};
            m.metNames  = {'a';'b'};
            m.metComps  = [1;1];
            m.metFormulas = {'H';'H'};
            m.metCharges  = [0;1];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.rxns      = {'R1'};
            m.rxnNames  = {'R1'};
            m.S         = sparse([-1;1]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0; m.b = zeros(2,1);
            m.grRules = {''}; m.genes = {}; m.rxnGeneMat = sparse(1,0);

            bs = getElementalBalance(m);
            testCase.verifyEqual(bs.chargeStatus(1), 0);
            testCase.verifyEqual(bs.chargeResidual(1), 1, 'AbsTol', 1e-9);
            % a -> b is elementally balanced (H on both sides)
            testCase.verifyEqual(bs.balanceStatus(1), 1);

            % Balancing the charge flips chargeStatus, not balanceStatus
            m.metCharges = [1;1];
            bs = getElementalBalance(m);
            testCase.verifyEqual(bs.chargeStatus(1), 1);
            testCase.verifyEqual(bs.chargeResidual(1), 0, 'AbsTol', 1e-9);
        end

        function getElementalBalanceChargeUnknownIsNotZero(testCase)
            % An unset charge on a participating metabolite makes the charge
            % balance unknown, and must not be reported as balanced. An unset
            % charge on a metabolite that does not participate must not leak
            % into the reaction's residual.
            m = testCase.model;
            m.mets      = {'a';'b';'spectator'};
            m.metNames  = {'a';'b';'spectator'};
            m.metComps  = [1;1;1];
            m.metFormulas = {'H';'H';'H'};
            m.metCharges  = [1;1;NaN];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.rxns      = {'R1'};
            m.rxnNames  = {'R1'};
            m.S         = sparse([-1;1;0]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0; m.b = zeros(3,1);
            m.grRules = {''}; m.genes = {}; m.rxnGeneMat = sparse(1,0);

            % The NaN belongs to a metabolite outside the reaction
            bs = getElementalBalance(m);
            testCase.verifyEqual(bs.chargeStatus(1), 1);
            testCase.verifyEqual(bs.chargeResidual(1), 0, 'AbsTol', 1e-9);

            % Now the NaN is on a participant: unknown, not balanced
            m.metCharges = [1;NaN;0];
            bs = getElementalBalance(m);
            testCase.verifyEqual(bs.chargeStatus(1), -1);
            testCase.verifyTrue(isnan(bs.chargeResidual(1)));
        end

        function getExchangeRxnsConsistent(testCase)
            [exch, idx] = getExchangeRxns(testCase.model);
            testCase.verifyClass(exch, 'cell');
            testCase.verifyNotEmpty(exch);
            testCase.verifyEqual(testCase.model.rxns(idx), exch);
        end

        function getExchangeRxnsWithBoundaryMets(testCase)
            % A model that still carries boundary metabolites (an
            % "unconstrained" field, as every model from importModel does)
            % takes a different branch than one where they have been
            % removed. Both have to return reaction indexes.
            m = testCase.taskTestModel();  % closeModel'd, so it has boundary mets
            testCase.assertTrue(isfield(m, 'unconstrained'));
            testCase.assertTrue(any(m.unconstrained ~= 0));

            [exch, idx] = getExchangeRxns(m);
            testCase.verifyEqual(sort(exch), {'R1'; 'R8'});
            testCase.verifyEqual(m.rxns(idx), exch);

            % R1 supplies a[s] ("=> a[s]"), R8 removes e[s] ("e[s] =>")
            testCase.verifyEqual(getExchangeRxns(m, 'in'), {'R1'});
            testCase.verifyEqual(getExchangeRxns(m, 'out'), {'R8'});
        end

        function getExchangeRxnsAgreeAcrossBoundaryRemoval(testCase)
            % Removing the boundary metabolites must not change which
            % reactions are exchange reactions.
            m = testCase.taskTestModel();
            withBoundary = getExchangeRxns(m);
            withoutBoundary = getExchangeRxns(simplifyModel(m));
            testCase.verifyEqual(sort(withoutBoundary), sort(withBoundary));
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
            % not be passed through as a logical array: this requires an
            % islogical check, since all() cannot distinguish the two cases.
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

        function parseFormulasRejectsPartialParse(testCase)
            % A formula that stops at an unrecognised element is not parsed,
            % even though the part before it was readable. Reporting it as
            % parsed would hand back a silently truncated composition.
            [elements, useMat, exitFlag] = parseFormulas({'C6H12O6';'C6Zz3'});
            testCase.verifyEqual(exitFlag, [1;-1]);
            testCase.verifyEqual(sum(useMat(2,:)), 0);
            % The good formula is unaffected
            cIdx = strcmp(elements.abbrevs, 'C');
            testCase.verifyEqual(useMat(1, cIdx), 6);
        end

        function parseFormulasSingleFormulaKeepsElements(testCase)
            % With one formula useMat is a row vector; the unused-element
            % pruning must still test each element separately.
            [elements, useMat] = parseFormulas({'C6H12O6'});
            testCase.verifySize(useMat, [1 numel(elements.abbrevs)]);
            testCase.verifyEqual(sort(elements.abbrevs(:)), {'C';'H';'O'});
        end

        function parseFormulasUnknownMassBlanksOnlyItsOwnFormula(testCase)
            % R and X have no mass. A formula that uses one gets MW = NaN;
            % formulas that do not must keep their weight.
            [~, ~, ~, MW] = parseFormulas({'C6H12O6';'C2R';'H2O';'CX'});
            testCase.verifyFalse(isnan(MW(1)));
            testCase.verifyTrue(isnan(MW(2)));
            testCase.verifyFalse(isnan(MW(3)));
            testCase.verifyTrue(isnan(MW(4)));
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

        function printFluxesRecognizesBoundaryMetExchanges(testCase)
            % A model closed via closeModel represents an exchange as
            % realMet <=> realMet[b] -- one reactant and one product --
            % which "no reactants or no products" alone does not recognise,
            % printing nothing for any of them.
            m = testCase.taskTestModel();
            [~,exchIdx] = getExchangeRxns(m,'all');
            testCase.assumeNotEmpty(exchIdx, 'Fixture has no exchange reactions after closeModel.');
            flux = ones(numel(m.rxns),1);
            out = evalc('printFluxes(m, flux, ''onlyExchange'', true);');
            testCase.verifySubstring(out, m.rxns{exchIdx(1)});
        end

        function getTransportRxnsExcludesBoundaryMetExchanges(testCase)
            % closeModel copies the real metabolite's own name onto its
            % boundary counterpart, so an exchange reaction has the same
            % "same name, different compartment" shape as a genuine
            % transport reaction.
            m = testCase.taskTestModel();
            [~,exchIdx] = getExchangeRxns(m,'all');
            testCase.assumeNotEmpty(exchIdx, 'Fixture has no exchange reactions after closeModel.');
            tr = getTransportRxns(m);
            testCase.verifyFalse(any(tr(exchIdx)));
        end

        function printFluxesRuns(testCase)
            testCase.assumeSolver('solveLP');
            sol = solveLP(testCase.model);
            out = evalc('printFluxes(testCase.model, sol.x, true)');
            testCase.verifyClass(out, 'char');
        end

        function printFluxesKeepsPercentInRxnName(testCase)
            % A rxnName containing a literal "%" must not be reinterpreted
            % as a format directive, which would truncate everything after
            % it instead of printing the flux line in full.
            m = testCase.model;
            m.rxnNames{1} = 'reaction with 50% yield';
            flux = zeros(numel(m.rxns),1);
            flux(1) = 1;
            out = evalc('printFluxes(m, flux, false)');
            testCase.verifySubstring(out, 'reaction with 50% yield');
        end

        function printModelStatsRuns(testCase)
            out = evalc('printModelStats(testCase.model)');
            testCase.verifyClass(out, 'char');
        end

        function printModelStatsKeepsPercentInModelName(testCase)
            % model.name/model.id/compNames/mets/rxns are spliced into
            % fprintf templates; a literal "%" in any of them must not be
            % reinterpreted as a format directive.
            m = testCase.model;
            m.name = 'test model with 50% coverage';
            m.S(1,:) = 0; % make mets{1} unused so it hits the errorText path
            m.metNames{1} = 'unused met 30% pure';
            out = evalc(['printModelStats(m, ''printModelIssues'', true, ' ...
                '''printDetails'', true);']);
            testCase.verifySubstring(out, 'test model with 50% coverage');
            testCase.verifySubstring(out, 'unused met 30% pure');
        end

    end
end
