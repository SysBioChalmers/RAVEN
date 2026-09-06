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

        function addGenesRavenPartialOverlapAddsRemainder(testCase)
            % When some genesToAdd.genes already exist in the model, the
            % others must still be added, and every parallel field must stay
            % aligned with the (trimmed) gene list.
            g.genes = {testCase.model.genes{1}; 'newgene1'};
            g.geneShortNames = {'existing','s1'};
            evalc('m2 = addGenesRaven(testCase.model, g);');
            testCase.verifyEqual(numel(m2.genes), numel(testCase.model.genes) + 1);
            testCase.verifyTrue(ismember('newgene1', m2.genes));
            testCase.verifyEqual(numel(m2.geneShortNames), numel(m2.genes));
        end

        function addMetsAddsMets(testCase)
            mta.metNames = {'newMetA','newMetB'};
            mta.compartments = {'c','e'};
            evalc('m2 = addMets(testCase.model, mta);');
            testCase.verifyEqual(numel(m2.mets), numel(testCase.model.mets) + 2);
        end

        function addRxnsAllowNewGenesKeepsGeneIdsIntact(testCase)
            % A gene id that merely contains "and"/"or" as a substring
            % (e.g. "band1") must not be shredded by the grRule parser used
            % to discover new genes.
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            r.grRules = 'band1 and orfeo2';
            evalc('m2 = addRxns(testCase.model, r, 2, ''c'', true, true);');
            testCase.verifyTrue(all(ismember({'band1','orfeo2'}, m2.genes)));
        end

        function addRxnsPrintsNewGeneIdWithPercent(testCase)
            % A new gene id containing "%" must not be truncated when
            % printed to the "New genes added" notice.
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            r.grRules = 'NEWGENE_50%_test';
            out = evalc('addRxns(testCase.model, r, 2, ''c'', true, true);');
            testCase.verifySubstring(out, 'NEWGENE_50%_test');
        end

        function addRxnsAddsRxn(testCase)
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            evalc('m2 = addRxns(testCase.model, r, 2, ''c'', true);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns) + 1);
            testCase.verifyTrue(ismember('newRxn1', m2.rxns));
        end

        function addRxnsKeepsSpontaneousAligned(testCase)
            % A model that already tracks spontaneous must keep it aligned
            % with rxns after adding reactions that don't specify it, and
            % accept an explicit spontaneous value for the new ones too.
            m = testCase.model;
            m.spontaneous = false(numel(m.rxns), 1);
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            evalc('m2 = addRxns(m, r, 2, ''c'', true);');
            testCase.verifyEqual(numel(m2.spontaneous), numel(m2.rxns));
            testCase.verifyFalse(m2.spontaneous(end));

            r2.rxns = 'newRxn2';
            r2.equations = '2-Oxoglutarate => TEST2';
            r2.spontaneous = true;
            evalc('m3 = addRxns(m, r2, 2, ''c'', true);');
            testCase.verifyEqual(numel(m3.spontaneous), numel(m3.rxns));
            testCase.verifyTrue(m3.spontaneous(end));
        end

        function addRxnsKeepsEquationsAligned(testCase)
            % A model that already tracks equations must keep it aligned
            % with rxns after adding a reaction, reusing the equation as
            % given rather than leaving the field short by one.
            m = testCase.model;
            m.equations = constructEquations(m);
            r.rxns = 'newRxn1';
            r.equations = '2-Oxoglutarate => TEST';
            evalc('m2 = addRxns(m, r, 2, ''c'', true);');
            testCase.verifyEqual(numel(m2.equations), numel(m2.rxns));
            testCase.verifySubstring(m2.equations{end}, 'TEST');
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

        function addRxnsGenesMetsKeepsPercentInAlreadyPresentNotice(testCase)
            % A rxn id containing "%" that is already present in the
            % target model must survive intact in the notice, not be
            % truncated by fprintf misreading it as a format directive.
            model = testCase.model;
            model.rxns{1} = 'RXN_50%_present';
            sourceModel = model;
            r.rxns = 'newRxn1';
            r.equations = [model.mets{1} ' => ' model.mets{2}];
            evalc('sourceModel = addRxns(sourceModel, r, 1, ''c'', true);');
            out = evalc(['addRxnsGenesMets(model, sourceModel, ' ...
                '{''RXN_50%_present'',''newRxn1''});']);
            testCase.verifySubstring(out, 'RXN_50%_present');
        end

        function addTransportAddsRxn(testCase)
            evalc(['m2 = addTransport(testCase.model, ''c'', ''e'', ' ...
                '{''6-phospho-D-glucono-1,5-lactone''}, false, false, ''tr_'');']);
            testCase.verifyGreaterThan(numel(m2.rxns), numel(testCase.model.rxns));
        end

        function addTransportAcceptsRowOrientedMetNames(testCase)
            % Regression test for SysBioChalmers/RAVEN#722: a metNames cell
            % array keeps whatever orientation the caller passed it in, and
            % a plain {'a','b'} literal is row-oriented. A single-metabolite
            % cell (as in addTransportAddsRxn above) cannot catch this ---
            % a 1x1 cell is the same shape whichever way you call it "row"
            % or "column" --- so this needs two or more metabolites, both
            % already present in the target compartment (the default
            % onlyToExisting=true path, which is where this crashed).
            metNames = {'CO2', 'Formate'}; %#ok<NASGU> % row-oriented on purpose
            evalc('m2 = addTransport(testCase.model, ''c'', {''e''}, ''metNames'', metNames);');
            testCase.verifyEqual(numel(m2.rxns), numel(testCase.model.rxns) + 2);
            testCase.verifyEqual(size(m2.rxnNames, 2), 1);
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

        function closeModelHandlesTwoColumnBAndMetNotes(testCase)
            % b(numel(b)+1)=0 is a linear index, which for an N-by-2 b
            % (net-production bounds) does not append a row; it errors.
            % metNotes must also be padded like the other optional
            % per-metabolite fields, or it goes out of sync with mets.
            m = tManipulation.twoMetModel();  % R1: a => b
            r.rxns = {'R2'}; r.equations = {'b =>'};  % sink, so closeModel adds a boundary met
            evalc('m = addRxns(m, r, 1, [], false);');
            m.b = [zeros(2,1) ones(2,1)];     % two-column b
            m.metNotes = {'note a';'note b'};
            m2 = closeModel(m);
            testCase.verifyEqual(size(m2.b,2), 2);
            testCase.verifyEqual(size(m2.b,1), numel(m2.mets));
            testCase.verifyEqual(numel(m2.metNotes), numel(m2.mets));
        end

        function closeModelDetectsScaledAndMultiMetSinks(testCase)
            % closeModel's boundary-reaction rule is "metabolites on only
            % one side" (matching getExchangeRxns), not "coefficients
            % summing to 1 in absolute value": a scaled single-metabolite
            % sink and a multi-metabolite one must both be detected and
            % closed, while a genuine two-sided reaction is left alone.
            m = tManipulation.twoMetModel();  % R1: a => b, a genuine reaction
            r.rxns = {'R2';'R3'};
            r.equations = {'2 a =>'; '0.5 a + 0.5 b =>'};
            evalc('m = addRxns(m, r, 1, [], false);');
            nMetsBefore = numel(m.mets);
            m2 = closeModel(m);

            % Exactly R2 and R3 are exchange-like, so exactly two boundary
            % metabolites are added, one per closed reaction.
            testCase.verifyEqual(numel(m2.mets) - nMetsBefore, 2);
            testCase.verifyEqual(numel(m2.metNames), numel(m2.mets));
            testCase.verifyEqual(numel(m2.metComps), numel(m2.mets));

            idx = getIndexes(m2, {'R1';'R2';'R3'}, 'rxns');
            boundaryComp = numel(m2.comps);
            touchesBoundary = full(any(m2.S(m2.metComps==boundaryComp, idx) ~= 0, 1));
            testCase.verifyEqual(touchesBoundary, [false true true]);
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

        function convertToIrrevReverseCopyInheritsSpontaneousAndPwys(testCase)
            % The _REV copy must carry the same spontaneous/pwys annotation
            % as the forward reaction, keeping both fields aligned with rxns.
            m = tManipulation.twoMetModel();
            m.rev = 1; m.lb = -500; m.ub = 1000;
            m.spontaneous = true;
            m.pwys = {'pathway1'};
            m2 = convertToIrrev(m);
            testCase.verifyEqual(numel(m2.spontaneous), numel(m2.rxns));
            testCase.verifyEqual(numel(m2.pwys), numel(m2.rxns));
            testCase.verifyTrue(m2.spontaneous(strcmp(m2.rxns,'R1_REV')));
            testCase.verifyEqual(m2.pwys{strcmp(m2.rxns,'R1_REV')}, 'pathway1');
        end

        function convertToIrrevRev2irrevPointsAtReverseCopy(testCase)
            % rev2irrev{origIdx}'s second element must be the actual
            % position of that reaction's reverse copy in irrevModel
            % (numOrigRxns+i), not just its rank among reversible
            % reactions.
            m = tManipulation.twoMetModel();       % R1: a -> b
            r.rxns = {'R2'}; r.equations = {'b <=> a'};
            evalc('m = addRxns(m, r, 1, [], false);');
            m.rev(2) = 1; m.lb(2) = -1000; m.ub(2) = 1000;
            [irrevModel, ~, rev2irrev] = convertToIrrev(m);
            pair = rev2irrev{2};
            testCase.verifyEqual(pair(1), 2);
            testCase.verifyEqual(irrevModel.rxns{pair(2)}, 'R2_REV');
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

        function mergeModelsKeepsSpontaneousAligned(testCase)
            % A model carrying spontaneous merged with one that doesn't
            % must keep the field aligned with rxns, defaulting the
            % other model's reactions to false rather than leaving it
            % short.
            a = tManipulation.namedMetModel('glc_c', 'A');
            a.spontaneous = true;
            b = tManipulation.namedMetModel('glucose_c', 'B');
            evalc('merged = mergeModels({a; b}, ''metParam'', ''mets'');');
            testCase.verifyEqual(numel(merged.spontaneous), numel(merged.rxns));
            testCase.verifyTrue(merged.spontaneous(strcmp(merged.rxns,'R_A')));
            testCase.verifyFalse(merged.spontaneous(strcmp(merged.rxns,'R_B')));
        end

        function copyToCompsDefaultCompOutsideAddsCompartment(testCase)
            % Adding a new compartment without specifying compOutside must
            % not error when the model already tracks compOutside.
            m = testCase.model;
            m.compOutside = repmat({''}, numel(m.comps), 1);
            evalc('m2 = copyToComps(m, {''p''}, ''ACKr'');');
            testCase.verifyEqual(numel(m2.compOutside), numel(m2.comps));
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

        function generateNewIdsEscapesRegexPrefix(testCase)
            % A prefix containing a regex metacharacter ('.') must be
            % matched literally: 'pX999' does not use the 'p.' prefix and
            % must not be picked up as if '.' were a wildcard.
            m.rxns = {'p.001'; 'pX999'};
            ids = generateNewIds(m, 'rxns', 'p.', 'quantity', 1);
            testCase.verifyEqual(ids{1}, 'p.002');
        end

        function mergeCompartmentsSingleComp(testCase)
            evalc('m2 = mergeCompartments(testCase.model);');
            testCase.verifyNumElements(m2.comps, 1);
        end

        function mergeCompartmentsWarnsWhenUnconstrainedMissing(testCase)
            % The warning's own text says it fires because there is no
            % unconstrained field to tell single-metabolite reactions apart
            % from real exchange reactions; the condition guarding it
            % checked the opposite, firing only when the field WAS present.
            m = testCase.model;
            testCase.verifyWarning(@() mergeCompartments(m,'deleteRxnsWithOneMet',true), ...
                'RAVEN:warning');
        end

        function mergeCompartmentsDropsStaleCompMiriams(testCase)
            % Merging collapses every compartment into one, so any
            % per-compartment MIRIAM annotation must not be left behind at the
            % original length (which would desync compMiriams from comps).
            m = testCase.model;
            m.compMiriams = cell(numel(m.comps),1);
            m.compMiriams{1}.name  = {'go'};
            m.compMiriams{1}.value = {'GO:0005737'};
            evalc('m2 = mergeCompartments(m);');
            if isfield(m2,'compMiriams')
                testCase.verifyNumElements(m2.compMiriams, numel(m2.comps));
            end
        end

        function mergeModelsReturnsStruct(testCase)
            modelB = removeReactions(testCase.model, testCase.model.rxns(1:5), true, true, true);
            evalc('merged = mergeModels({testCase.model; modelB});');
            testCase.verifyClass(merged, 'struct');
            testCase.verifyNotEmpty(merged.rxns);
        end

        function mergeModelsSelfCollisionGetsUniqueIds(testCase)
            % A source model with two internally duplicate reaction ids
            % ('R1' twice) that both collide with the growing merged model
            % must each be renamed to a distinct id, not the same one.
            a = tManipulation.twoMetModel();
            b = tManipulation.duplicateRxnModel('M2');
            evalc('merged = mergeModels({a; b});');
            testCase.verifyEqual(numel(unique(merged.rxns)), numel(merged.rxns));
            newIds = merged.rxns(2:3);
            testCase.verifyTrue(all(startsWith(newIds, 'R1_M2')));
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

        function removeBadRxnsSeedIsReproducible(testCase)
            % Which reaction is removed among several equally-valid
            % candidates is randomised; a given seed must make the choice
            % (and thus the result) reproducible across runs.
            evalc(['[~, r1] = removeBadRxns(testCase.model, ''rxnRules'', 3, ''seed'', 42);' ...
                '[~, r2] = removeBadRxns(testCase.model, ''rxnRules'', 3, ''seed'', 42);']);
            testCase.verifyEqual(r1, r2);
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

        function replaceMetsVerboseKeepsPercentInRxnId(testCase)
            % A reaction id containing "%" that is reported by 'verbose'
            % must survive intact, not be truncated by fprintf misreading
            % it as a format directive.
            m = struct();
            m.id='t'; m.rxns={'R_50%_test';'R2'}; m.rxnNames=m.rxns;
            m.mets={'x';'y';'z'}; m.metNames=m.mets; m.metComps=[1;1;1];
            m.comps={'c'}; m.compNames={'c'};
            m.S=sparse([-1 0; 0 -1; 1 1]); % R_50%_test: x=>z   R2: y=>z
            m.lb=[0;0]; m.ub=[1000;1000]; m.rev=[0;0]; m.c=[0;0]; m.b=zeros(3,1);
            m.genes={}; m.grRules={'';''}; m.rxnGeneMat=sparse(2,0);
            out = evalc('replaceMets(m, ''x'', ''y'', ''verbose'', true);');
            testCase.verifySubstring(out, 'R_50%_test');
        end

        function replaceMetsByIdAddsRatherThanOverwrites(testCase)
            % A reaction where the replacement metabolite is already itself
            % a participant must keep that contribution: x+y=> must become
            % 2y=>, not just y=>, once x is replaced by y.
            m = struct();
            m.id='t'; m.rxns={'R1';'R2'}; m.rxnNames=m.rxns;
            m.mets={'x';'y';'z'}; m.metNames=m.mets; m.metComps=[1;1;1];
            m.comps={'c'}; m.compNames={'c'};
            m.S=sparse([0 -1; 0 -1; -1 0]); % R1: z=>   R2: x+y=>
            m.lb=[0;0]; m.ub=[1000;1000]; m.rev=[0;0]; m.c=[0;0]; m.b=zeros(3,1);
            m.genes={}; m.grRules={'';''}; m.rxnGeneMat=sparse(2,0);
            evalc('m2 = replaceMets(m, ''x'', ''y'', ''identifiers'', true);');
            yRow = strcmp(m2.mets,'y');
            r2 = strcmp(m2.rxns,'R2');
            testCase.verifyEqual(full(m2.S(yRow,r2)), -2);
        end

        function replaceMetsByNameKeepsBShapeAndRuns(testCase)
            % A model with a two-column b (net-production bounds) must
            % keep that shape after the metabolites-with-duplicate-name
            % merge, and the final contractModel call must not error from
            % a stale post-deletion metabolite index or a disabled
            % distReverse.
            m = struct();
            m.id='t'; m.rxns={'R1';'R2'}; m.rxnNames=m.rxns;
            m.mets={'ox1';'ox2';'w'}; m.metNames={'oxygen';'o2';'w'};
            m.metComps=[1;1;1]; m.comps={'c'}; m.compNames={'c'};
            m.S=sparse([-1 0; 0 -1; 1 1]); % R1: oxygen=>w   R2: o2=>w
            m.lb=[0;0]; m.ub=[1000;1000]; m.rev=[0;0]; m.c=[0;0];
            m.b=[zeros(3,1) ones(3,1)]; % two-column b
            m.genes={}; m.grRules={'';''}; m.rxnGeneMat=sparse(2,0);
            evalc('m2 = replaceMets(m, ''oxygen'', ''o2'');');
            testCase.verifyEqual(size(m2.b,2), 2);
            testCase.verifyEqual(size(m2.b,1), numel(m2.mets));
            testCase.verifyEqual(numel(m2.mets), 2); % oxygen and o2 merged
        end

        function setExchangeBoundsFindsAllMultiExchangeMets(testCase)
            % Every metabolite exchanged by more than one reaction must be
            % reported, not just whichever one happens to line up between
            % two differently-sized index ranges compared directly against
            % each other.
            m = struct();
            m.id='t';
            m.rxns={'r1';'r2';'r3';'r4';'r5';'r6';'r7'}; m.rxnNames=m.rxns;
            m.mets={'met1';'met2';'met3';'met4';'met5'};
            m.metNames={'MetOne';'MetTwo';'MetThree';'MetFour';'MetFive'};
            m.metComps=[1;1;1;1;1]; m.comps={'c'}; m.compNames={'c'};
            % r1:=>met5  r2:=>met1  r3:=>met2  r4:met5=>  r5:=>met3  r6:=>met4  r7:met3=>
            m.S = sparse(5,7);
            m.S(5,1)=1; m.S(1,2)=1; m.S(2,3)=1; m.S(5,4)=-1; m.S(3,5)=1; m.S(4,6)=1; m.S(3,7)=-1;
            m.lb=-1000*ones(7,1); m.ub=1000*ones(7,1); m.rev=ones(7,1); m.c=zeros(7,1); m.b=zeros(5,1);
            m.genes={}; m.grRules=repmat({''},7,1); m.rxnGeneMat=sparse(7,0);
            txt = evalc('setExchangeBounds(m);');
            testCase.verifySubstring(txt, 'MetThree');
            testCase.verifySubstring(txt, 'MetFive');
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

        function simplifyModelKeepsMetsConstrainedByB(testCase)
            % R2 is constrained to zero flux, so deleteZeroInterval removes
            % it and leaves C in no reaction at all. A non-zero b is a
            % boundary condition on C, and removing the metabolite takes
            % that row of b with it, so whatever solves the reduced model
            % is no longer asked to produce C.
            m = struct();
            m.id = 'test';
            m.rxns = {'R1';'R2'}; m.rxnNames = m.rxns;
            m.mets = {'A';'B';'C'}; m.metNames = m.mets; m.metComps = [1;1;1];
            m.comps = {'c'}; m.compNames = m.comps;
            m.S = sparse([-1 0; 1 -1; 0 1]);   % R1: A => B, R2: B => C
            m.lb = [0;0]; m.ub = [1000;0]; m.rev = [0;0]; m.c = [0;0];
            m.genes = {}; m.grRules = {'';''}; m.rxnGeneMat = sparse(2,0);
            m.b = zeros(3,2);
            m.b(3,:) = [1 1];                  % C must be produced, 1 unit

            evalc('[reduced,~,deletedMets] = simplifyModel(m,false,false,true);');
            testCase.verifyTrue(ismember('C', reduced.mets));
            testCase.verifyFalse(ismember('C', deletedMets));
            testCase.verifyEqual(reduced.b(strcmp(reduced.mets,'C'),:), [1 1]);

            % Without that constraint C is dead weight and is removed.
            m.b(3,:) = [0 0];
            evalc('reduced = simplifyModel(m,false,false,true);');
            testCase.verifyFalse(ismember('C', reduced.mets));
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

        function sortReactionOrderSeedIsReproducible(testCase)
            % The local search proposes swaps randomly; a given seed must
            % make its result reproducible across runs.
            m = testCase.model;
            m.subSystems = repmat({{'ALL'}}, numel(m.rxns), 1);
            m1 = sortModel(m, 'sortReversible', false, 'sortReactionOrder', true, 'seed', 5);
            m2 = sortModel(m, 'sortReversible', false, 'sortReactionOrder', true, 'seed', 5);
            testCase.verifyEqual(m1.rxns, m2.rxns);
        end

        function sortReactionOrderUsesSubsystemsOwnColumns(testCase)
            % sortReactionOrder must score and reorder a subsystem's own
            % reactions, not whichever columns happen to occupy the first
            % nRxns positions of the whole model: with a chain A->B->C->D
            % split across R2 (A=>B), R3 (B=>C) and R1 (C=>D) -- placed
            % after two unrelated filler reactions so the subsystem is NOT
            % at the start of model.rxns -- the only production-before-
            % consumption order is R2, then R3, then R1.
            m = struct();
            m.id='t'; m.rxns={'FILLER1';'FILLER2';'R1';'R2';'R3'}; m.rxnNames=m.rxns;
            m.mets={'fa';'fb';'fc';'fd';'A';'B';'C';'D'}; m.metNames=m.mets;
            m.metComps=ones(8,1); m.comps={'c'}; m.compNames={'c'};
            S=zeros(8,5);
            S(1,1)=-1; S(2,1)=1;  % FILLER1: fa=>fb
            S(3,2)=-1; S(4,2)=1;  % FILLER2: fc=>fd
            S(7,3)=-1; S(8,3)=1;  % R1: C=>D
            S(5,4)=-1; S(6,4)=1;  % R2: A=>B
            S(6,5)=-1; S(7,5)=1;  % R3: B=>C
            m.S=sparse(S);
            m.lb=zeros(5,1); m.ub=ones(5,1)*1000; m.rev=zeros(5,1); m.c=zeros(5,1);
            m.b=zeros(8,1);
            m.genes={}; m.grRules=repmat({''},5,1); m.rxnGeneMat=sparse(5,0);
            m.subSystems={{};{};{'CHAIN'};{'CHAIN'};{'CHAIN'}};

            rng(1);
            m2=sortModel(m,'sortReversible',false,'sortReactionOrder',true);
            posR1=find(strcmp(m2.rxns,'R1'));
            posR2=find(strcmp(m2.rxns,'R2'));
            posR3=find(strcmp(m2.rxns,'R3'));
            testCase.verifyLessThan(posR2, posR3);
            testCase.verifyLessThan(posR3, posR1);
        end

        function standardizeGrRulesReturnsRules(testCase)
            evalc('grRules = standardizeGrRules(testCase.model);');
            testCase.verifyNumElements(grRules, numel(testCase.model.rxns));
        end

        function findPotentialErrorsFlagsOnlyNonDnf(testCase)
            m = struct();
            m.rxns = {'R1';'R2';'R3';'R4'};
            m.grRules = {'((G1 and G2) or G3)'      % DNF, just bracketed
                         '(G1 or G2) or (G3 and G4)' % DNF
                         'G1 or G2'                  % DNF
                         '(G1 or G2) and (G3 or G4)'};% genuinely non-DNF
            issues = findPotentialErrors(m);
            testCase.verifyEqual(vertcat(issues.index), 4);
        end

        function findPotentialErrorsReportsUnparseable(testCase)
            m = struct();
            m.rxns = {'R1'};
            m.grRules = {'(G1 and G2'};
            issues = findPotentialErrors(m);
            testCase.verifyNumElements(issues, 1);
            testCase.verifySubstring(issues(1).reason, 'Cannot be parsed');
        end

        function standardizeGrRulesRepairsBracketedDnf(testCase)
            % A rule that is DNF but redundantly bracketed must be repaired,
            % not skipped: standardizeGrRules leaves flagged rules alone, so a
            % false positive from findPotentialErrors silently prevents repair.
            m = struct();
            m.rxns = {'R1'};
            m.grRules = {'((G1 and G2) or G3)'};
            m.genes = {'G1';'G2';'G3'};
            m.rxnGeneMat = sparse([1 1 1]);
            [grRules,~,indexes2check] = standardizeGrRules(m, true);
            testCase.verifyEmpty(indexes2check);
            testCase.verifyEqual(grRules{1}, '(G1 and G2) or G3');
        end

        function standardizeGrRulesKeepsPercentInWarning(testCase)
            % A rxn id containing a literal "%" must survive intact in the
            % "potentially problematic relationships" warning, not be
            % truncated by sprintf/warning misreading it as a directive.
            m.rxns = {'RXN_50%_TEST'};
            m.grRules = {'(G1 or G2) and G3'};
            m.genes = {'G1';'G2';'G3'};
            m.rxnGeneMat = sparse([1 1 1]);
            lastwarn('');
            evalc('standardizeGrRules(m);');
            msg = lastwarn();
            testCase.verifySubstring(msg, 'RXN_50%_TEST');
        end

        function removeGenesMatchesWholeGeneIds(testCase)
            % Gene "10" is a prefix of "100". Removing it must not take "100"
            % with it, which an unanchored substring search would.
            m = testCase.gprTestModel('10 or 100', {'10';'100'}, [1 1]);
            r = removeGenes(m, {'10'});
            testCase.verifyEqual(r.grRules{1}, '100');
            testCase.verifyEqual(r.ub(1), 1000);
        end

        function removeGenesDropsWholeComplex(testCase)
            % A complex missing a subunit cannot form; the other isozyme lives.
            m = testCase.gprTestModel('(G1 and G2) or G3', {'G1';'G2';'G3'}, [1 1 1]);
            r = removeGenes(m, {'G1'});
            testCase.verifyEqual(r.grRules{1}, 'G3');
        end

        function removeGenesBlocksWhenNoEnzymeLeft(testCase)
            m = testCase.gprTestModel('G1 and G2', {'G1';'G2'}, [1 1]);
            r = removeGenes(m, {'G1'});
            testCase.verifyEmpty(r.grRules{1});
            testCase.verifyEqual(r.lb(1), 0);
            testCase.verifyEqual(r.ub(1), 0);
        end

        function expandModelKeepsMandatorySubunit(testCase)
            % "g1 and (g2 or g3)" is two isozymes, both needing g1. Stripping
            % brackets and splitting on ' or ' loses g1 from the second.
            m = testCase.gprTestModel('g1 and (g2 or g3)', {'g1';'g2';'g3'}, [1 1 1]);
            e = expandModel(m);
            testCase.verifyEqual(sort(e.grRules), {'g1 and g2';'g1 and g3'});
        end

        function expandModelCopiesSpontaneousAndPwys(testCase)
            % Each isozyme copy created by splitting an OR rule must inherit
            % the source reaction's spontaneous/pwys annotation, keeping both
            % fields aligned with rxns.
            m = testCase.gprTestModel('g1 or g2', {'g1';'g2'}, [1 1]);
            m.spontaneous = true;
            m.pwys = {'pathway1'};
            e = expandModel(m);
            testCase.verifyEqual(numel(e.spontaneous), numel(e.rxns));
            testCase.verifyEqual(numel(e.pwys), numel(e.rxns));
            testCase.verifyTrue(all(e.spontaneous));
            testCase.verifyTrue(all(strcmp(e.pwys, 'pathway1')));
        end

        function expandModelCopiesRxnScores(testCase)
            % Each isozyme copy created by splitting an OR rule must
            % inherit the source reaction's rxnScores, keeping it aligned
            % with rxns.
            m = testCase.gprTestModel('g1 or g2', {'g1';'g2'}, [1 1]);
            m.rxnScores = 2.5;
            e = expandModel(m);
            testCase.verifyEqual(numel(e.rxnScores), numel(e.rxns));
            testCase.verifyTrue(all(e.rxnScores == 2.5));
        end

        function expandModelDistributesBothSides(testCase)
            % Two or:s, but four isozymes.
            m = testCase.gprTestModel('(g1 or g2) and (g3 or g4)', ...
                {'g1';'g2';'g3';'g4'}, [1 1 1 1]);
            [e, rxnToCheck] = expandModel(m);
            testCase.verifyEqual(sort(e.grRules), ...
                {'g1 and g3';'g1 and g4';'g2 and g3';'g2 and g4'});
            testCase.verifyEqual(rxnToCheck, {'R1'});
        end

        function expandModelLeavesDnfAlone(testCase)
            % An OR of complexes expands without needing distributivity, so it
            % must not be reported as needing a check.
            m = testCase.gprTestModel('(g1 and g2) or (g3 and g4)', ...
                {'g1';'g2';'g3';'g4'}, [1 1 1 1]);
            [e, rxnToCheck] = expandModel(m);
            testCase.verifyEqual(sort(e.grRules), {'g1 and g2';'g3 and g4'});
            testCase.verifyEmpty(rxnToCheck);
        end

    end

    methods (Access = private)
        function m = gprTestModel(~, grRule, genes, rxnGeneRow)
            % Smallest model that removeGenes/expandModel will operate on.
            m = struct();
            m.rxns = {'R1'}; m.rxnNames = {'R1'};
            m.mets = {'A';'B'}; m.metNames = {'A';'B'}; m.metComps = [1;1];
            m.comps = {'c'}; m.compNames = {'c'};
            m.S = sparse([-1;1]); m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0;
            m.b = [0;0];
            m.genes = genes;
            m.grRules = {grRule};
            m.rxnGeneMat = sparse(rxnGeneRow);
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

        function m = duplicateRxnModel(modelId)
            % Two reactions, both named 'R1' internally. mergeModels'
            % own docstring says duplicate reaction ids "might appear" in
            % a model, so this is a real (if unusual) input shape to
            % guard the merge against, not a contrived one.
            m = struct();
            m.id        = modelId;
            m.rxns      = {'R1'; 'R1'};
            m.rxnNames  = {'R1a'; 'R1b'};
            m.mets      = {'x'; 'y'; 'z'};
            m.metNames  = m.mets;
            m.metComps  = [1;1;1];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-1 0; 1 -1; 0 1]);
            m.lb = [0;0]; m.ub = [1000;1000]; m.rev = [0;0]; m.c = [0;0]; m.b = zeros(3,1);
            m.genes = {}; m.grRules = {'';''}; m.rxnGeneMat = sparse(2,0);
            m.metFormulas = {'C';'C';'C'};
        end

    end
end
