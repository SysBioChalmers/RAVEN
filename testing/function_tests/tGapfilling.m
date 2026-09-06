classdef tGapfilling < RavenTestCase
% tGapfilling  Tests for the gap-analysis and gap-filling functions in gapfilling/.

    methods (Test)

        function canExchangeConsumeReturnsLogical(testCase)
            out = canExchange(testCase.model, 'consume', testCase.model.mets(1:3));
            testCase.verifyClass(out, 'logical');
            testCase.verifyNumElements(out, 3);
        end

        function canExchangeProduceReturnsLogical(testCase)
            out = canExchange(testCase.model, 'produce', testCase.model.mets(1:3));
            testCase.verifyClass(out, 'logical');
            testCase.verifyNumElements(out, 3);
        end

        function canExchangeInvalidDirectionErrors(testCase)
            testCase.verifyError( ...
                @() canExchange(testCase.model, 'neither'), ...
                'RAVEN:badInput');
        end

        function checkProductionReturnsIndices(testCase)
            evalc('notProduced = checkProduction(testCase.model);');
            testCase.verifyClass(notProduced, 'double');
        end

        function findLeakMetaboliteProduceReturnsSolution(testCase)
            evalc('[sol, metabolite] = findLeakMetabolite(testCase.model, ''produce'');');
            testCase.verifyNotEmpty(sol);
        end

        function findLeakMetaboliteConsumeReturnsSolution(testCase)
            evalc('[sol, metabolite] = findLeakMetabolite(testCase.model, ''consume'');');
            testCase.verifyNotEmpty(sol);
        end

        function findLeakMetaboliteInvalidDirectionErrors(testCase)
            testCase.verifyError( ...
                @() findLeakMetabolite(testCase.model, 'neither'), ...
                'RAVEN:badInput');
        end

        function deprecatedMakeSomethingStillWorks(testCase)
            % The deprecated/ wrapper must keep returning what it always did.
            % Its warning is once-per-session, so it is not asserted here
            % (ordering across tests would decide whether it fires); the
            % warning itself is covered by tUtils/deprecationWarning*.
            evalc('[sol, metabolite] = makeSomething(testCase.model);');
            testCase.verifyNotEmpty(sol);
        end

        function deprecatedConsumeSomethingStillWorks(testCase)
            evalc('[sol, metabolite] = consumeSomething(testCase.model);');
            testCase.verifyNotEmpty(sol);
        end

        function deprecatedConsumeSomethingKeepsItsOwnArgumentOrder(testCase)
            % consumeSomething's positional order has no allowExcretion,
            % unlike findLeakMetabolite's. Passing its 5th argument
            % (ignoreIntBounds) must reach ignoreIntBounds, not params --
            % forwarding varargin verbatim used to shift it silently.
            evalc(['[sol, metabolite] = consumeSomething(testCase.model, ' ...
                '[], false, false, [], true);']);
            testCase.verifyNotEmpty(sol);
        end

        function fillGapsProducesModel(testCase)
            testCase.assumeMILPSolver();
            modelDB = testCase.model; modelDB.id = 'DB';   % template
            gapModel = removeReactions(modelDB, (1:10));
            gapModel.id = 'gapModel';                      % must differ from template
            evalc('[~,~,~,newModel] = fillGaps(gapModel, modelDB);');
            testCase.verifyClass(newModel, 'struct');
        end

        function gapReportRuns(testCase)
            evalc('noFluxRxns = gapReport(testCase.model);');
            testCase.verifyClass(noFluxRxns, 'cell');
        end

        function gapReportKeepsPercentInModelName(testCase)
            % A model.name containing "%" must survive intact in the
            % report header, not be truncated by fprintf misreading it as
            % a format directive.
            m = testCase.model;
            m.name = 'ecoli 50% subset';
            out = evalc('gapReport(m);');
            testCase.verifySubstring(out, 'ecoli 50% subset');
        end

        function fitTasksProducesModel(testCase)
            testCase.assumeMILPSolver();
            refModel = testCase.taskTestModel(); refModel.id = 'DB';
            gapModel = removeReactions(refModel, {'R2'}); gapModel.id = 'testModel';
            task = testCase.taskTestStruct();
            evalc('[outModel, addedRxns] = fitTasks(gapModel, refModel, [], true, [], task);');
            testCase.verifyClass(outModel, 'struct');
            testCase.verifyTrue(ismember('R2', outModel.rxns));
        end

        function gapFillFastCoreReturnsLogical(testCase)
            % gapFillFastCore should return a logical vector the same length as model.rxns.
            model = testCase.model;
            coreIdx = 1;   % any single reaction as the core
            active = gapFillFastCore(model, coreIdx, 1e-4);
            testCase.verifyClass(active, 'logical');
            testCase.verifyNumElements(active, numel(model.rxns));
            testCase.verifyTrue(active(coreIdx));  % core reaction must be active
        end

        function gapFillFastCoreRejectsTriviallySelfCancelingReversibleCore(testCase)
            % A reversible core reaction touching only its own,
            % otherwise-unused metabolites cannot carry any real
            % steady-state flux (mass balance on those metabolites forces
            % v=0). Forcing both its forward and reverse irreversible
            % copies to >= epsilon must not let them satisfy the core
            % requirement by canceling each other out.
            model = testCase.model;
            r.rxns = {'isolatedRev'};
            r.equations = {'newA[c] <=> newB[c]'};
            evalc('model = addRxns(model, r, 3, ''c'', true);');
            coreIdx = getIndexes(model,'isolatedRev','rxns');
            active = gapFillFastCore(model, coreIdx, 1e-4);
            testCase.verifyFalse(active(coreIdx));
        end

        function gapFillSwiftCoreReturnsLogical(testCase)
            % gapFillSwiftCore returns same shape as gapFillFastCore.
            model = testCase.model;
            coreIdx = 1;
            active = gapFillSwiftCore(model, coreIdx, 1e-4);
            testCase.verifyClass(active, 'logical');
            testCase.verifyNumElements(active, numel(model.rxns));
            testCase.verifyTrue(active(coreIdx));
        end

        function gapFillSwiftCoreRejectsTriviallySelfCancelingReversibleCore(testCase)
            % Same issue as gapFillFastCore's identically-named test: a
            % reversible core reaction touching only its own,
            % otherwise-unused metabolites cannot carry any real
            % steady-state flux, so forcing both its forward and reverse
            % irreversible copies to >= epsilon must not let them satisfy
            % the core requirement by canceling each other out.
            model = testCase.model;
            r.rxns = {'isolatedRev'};
            r.equations = {'newA[c] <=> newB[c]'};
            evalc('model = addRxns(model, r, 3, ''c'', true);');
            coreIdx = getIndexes(model,'isolatedRev','rxns');
            active = gapFillSwiftCore(model, coreIdx, 1e-4);
            testCase.verifyFalse(active(coreIdx));
        end

        function fillGapsIdentifiesOwnRxnsWhenRxnFromPreset(testCase)
            % model.rxnFrom, as e.g. getModelFromHomology output already
            % carries it, must not stop fillGaps from recognising the
            % model's own reactions when checking which of them regain
            % flux via the template.
            testCase.assumeMILPSolver();
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, (1:10));
            gapModel.id = 'gapModel';
            gapModel.rxnFrom = repmat({'someTemplate'}, numel(gapModel.rxns), 1);
            evalc('[newConnected,~,~,newModel] = fillGaps(gapModel, modelDB);');
            testCase.verifyNotEmpty(newConnected);
            testCase.verifyTrue(all(ismember(gapModel.rxns, newModel.rxns)));
        end

        function gapFillFastLPReturnsAddedRxns(testCase)
            % gapFillFastLP should identify database reactions that rescue blocked draft reactions.
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('[addedRxns, newModel, cannotConnect] = gapFillFastLP(gapModel, modelDB, ''verbose'', false);');
            testCase.verifyClass(addedRxns, 'cell');
            testCase.verifyClass(newModel, 'struct');
            testCase.verifyClass(cannotConnect, 'cell');
            % Added reactions must all be present in the returned model and be
            % universal-derived (not original draft reactions). mergeModels
            % appends the template model id to deduplicate reaction IDs shared by
            % draft and DB, so added IDs may carry a '_DB' suffix rather than
            % matching modelDB.rxns verbatim.
            if ~isempty(addedRxns)
                testCase.verifyTrue(all(ismember(addedRxns, newModel.rxns)));
                testCase.verifyFalse(any(ismember(addedRxns, gapModel.rxns)));
            end
        end

        function gapFillSwiftLPReturnsAddedRxns(testCase)
            % swiftLP variant should produce consistent results with fastLP.
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('[addedRxns, newModel, ~] = gapFillFastLP(gapModel, modelDB, ''variant'', ''swift'', ''verbose'', false);');
            testCase.verifyClass(addedRxns, 'cell');
            testCase.verifyClass(newModel, 'struct');
        end

        function gapFillTopologicalIdentifiesGaps(testCase)
            % gapFillTopological should return a struct with the expected fields.
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('result = gapFillTopological(gapModel, modelDB, ''verbose'', false);');
            testCase.verifyClass(result, 'struct');
            testCase.verifyTrue(isfield(result, 'reachableMets'));
            testCase.verifyTrue(isfield(result, 'blockedMets'));
            testCase.verifyTrue(isfield(result, 'candidateRxns'));
            testCase.verifyTrue(isfield(result, 'pruningFraction'));
            testCase.verifyClass(result.reachableMets, 'logical');
            testCase.verifyNumElements(result.reachableMets, numel(gapModel.mets));
        end

        function gapFillTopologicalReachesThroughInExchange(testCase)
            % A -> B -> C fed by an 'in' exchange on A. addExchangeRxns writes
            % 'in' as "=> A", which has lb = 0, so a seed test of lb < 0 finds
            % no uptake at all and every metabolite comes back blocked.
            m = testCase.chainModel();
            m = addExchangeRxns(m, 'in', {'A'});
            evalc(['result = gapFillTopological(m, m, ''targets'', ' ...
                '{''A'',''B'',''C''}, ''verbose'', false);']);
            testCase.verifyTrue(all(result.reachableMets));
            testCase.verifyEmpty(result.blockedMets);
        end

        function gapFillTopologicalReachesThroughReversibleExchange(testCase)
            % "A <=>" supplies A too, but is 'reverse' rather than 'uptake'.
            m = testCase.chainModel();
            m = addExchangeRxns(m, 'both', {'A'});
            evalc(['result = gapFillTopological(m, m, ''targets'', ' ...
                '{''A'',''B'',''C''}, ''verbose'', false);']);
            testCase.verifyTrue(all(result.reachableMets));
        end

        function gapFillTopologicalBlocksWithoutExchange(testCase)
            % Without any exchange nothing is producible: the positive control
            % for the two tests above.
            m = testCase.chainModel();
            evalc(['result = gapFillTopological(m, m, ''seeds'', {}, ''targets'', ' ...
                '{''A'',''B'',''C''}, ''verbose'', false);']);
            testCase.verifyFalse(any(result.reachableMets));
            testCase.verifyNumElements(result.blockedMets, 3);
        end

        function gapFillMILPRepairsGrowth(testCase)
            testCase.assumeMILPSolver();
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('[addedRxns, reversedRxns, newModel, exitFlag] = gapFillMILP(gapModel, modelDB, ''verbose'', false);');
            testCase.verifyClass(addedRxns, 'cell');
            testCase.verifyClass(reversedRxns, 'cell');
            testCase.verifyClass(newModel, 'struct');
            testCase.verifyEqual(exitFlag, 1);
            % The repaired model should be able to produce objective flux.
            sol = solveLP(newModel);
            testCase.verifyNotEmpty(sol.f);
            testCase.verifyGreaterThan(sol.f, 0);
        end

        function gapFillMILPHandlesUnboundedUniversalReaction(testCase)
            % A universal reaction with ub=Inf (RAVEN's own default when no
            % upper-bound annotation is set) must not put an Inf
            % coefficient into the coupling constraint matrix -- SCIP
            % rejects that outright, and Inf*0 is NaN even where it
            % wouldn't be rejected.
            testCase.assumeMILPSolver();
            gapModel = struct();
            gapModel.id='gapModel'; gapModel.rxns={'R1';'ExB'}; gapModel.rxnNames=gapModel.rxns;
            gapModel.mets={'a';'b'}; gapModel.metNames=gapModel.mets; gapModel.metComps=[1;1];
            gapModel.comps={'c'}; gapModel.compNames={'c'};
            gapModel.S=sparse([-1 0; 1 -1]); % R1: a=>b   ExB: b=>
            gapModel.lb=[0;0]; gapModel.ub=[1000;1000]; gapModel.rev=[0;0];
            gapModel.c=[1;0]; gapModel.b=zeros(2,1);
            gapModel.genes={}; gapModel.grRules={'';''}; gapModel.rxnGeneMat=sparse(2,0);

            modelDB = struct();
            modelDB.id='DB'; modelDB.rxns={'ExA'}; modelDB.rxnNames=modelDB.rxns;
            modelDB.mets={'a'}; modelDB.metNames=modelDB.mets; modelDB.metComps=1;
            modelDB.comps={'c'}; modelDB.compNames={'c'};
            modelDB.S=sparse(1,1); modelDB.S(1,1)=1;
            modelDB.lb=0; modelDB.ub=Inf; modelDB.rev=0; modelDB.c=0; modelDB.b=0;
            modelDB.genes={}; modelDB.grRules={''}; modelDB.rxnGeneMat=sparse(1,0);

            evalc(['[addedRxns,~,newModel,exitFlag] = gapFillMILP(gapModel, modelDB, ' ...
                '''verbose'', false);']);
            testCase.verifyEqual(exitFlag, 1);
            testCase.verifyEqual(addedRxns, {'ExA'});
            sol = solveLP(newModel);
            testCase.verifyGreaterThan(sol.f, 0);
        end

        function gapFillMILPKeepsRenamedUniversalReaction(testCase)
            % A universal reaction whose id collides with an unrelated
            % draft reaction gets renamed by mergeModels; that renamed
            % reaction must still end up in addedRxns and in newModel, not
            % be silently dropped by a lookup keyed on its original name.
            testCase.assumeMILPSolver();
            gapModel = struct();
            gapModel.id='gapModel'; gapModel.rxns={'R1';'Rname'}; gapModel.rxnNames=gapModel.rxns;
            gapModel.mets={'a';'b'}; gapModel.metNames=gapModel.mets; gapModel.metComps=[1;1];
            gapModel.comps={'c'}; gapModel.compNames={'c'};
            gapModel.S=sparse([-1 0; 1 -1]); % R1: a=>b   Rname (draft): b=>, irrelevant to growth
            gapModel.lb=[0;0]; gapModel.ub=[1000;1000]; gapModel.rev=[0;0];
            gapModel.c=[1;0]; gapModel.b=zeros(2,1);
            gapModel.genes={}; gapModel.grRules={'';''}; gapModel.rxnGeneMat=sparse(2,0);

            modelDB = struct();
            % Same id "Rname" as the draft's, but a completely different
            % reaction: this is the one actually needed to enable growth.
            modelDB.id='DB'; modelDB.rxns={'Rname'}; modelDB.rxnNames=modelDB.rxns;
            modelDB.mets={'a'}; modelDB.metNames=modelDB.mets; modelDB.metComps=1;
            modelDB.comps={'c'}; modelDB.compNames={'c'};
            modelDB.S=sparse(1,1); modelDB.S(1,1)=1; % =>a
            modelDB.lb=0; modelDB.ub=1000; modelDB.rev=0; modelDB.c=0; modelDB.b=0;
            modelDB.genes={}; modelDB.grRules={''}; modelDB.rxnGeneMat=sparse(1,0);

            evalc(['[addedRxns,~,newModel,exitFlag] = gapFillMILP(gapModel, modelDB, ' ...
                '''verbose'', false);']);
            testCase.verifyEqual(exitFlag, 1);
            testCase.verifyNotEmpty(addedRxns);
            sol = solveLP(newModel);
            testCase.verifyGreaterThan(sol.f, 0);
        end

        function gapFillMILPReversesDirectionality(testCase)
            % gapFillMILP should reverse a reaction whose directionality is wrong.
            testCase.assumeMILPSolver();
            model = testCase.model;
            modelDB = model; modelDB.id = 'DB';
            gapModel = removeReactions(model, model.rxns(1:3));
            gapModel.id = 'gapModel';
            evalc('[~, reversedRxns, ~, ~] = gapFillMILP(gapModel, modelDB, ''verbose'', false);');
            testCase.verifyClass(reversedRxns, 'cell');
        end

        function fillGapsDispatchesFastLP(testCase)
            % fillGaps with 'algorithm','fastLP' should reach gapFillFastLP.
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('[~,cannotConnect,addedRxns,newModel,exitFlag] = fillGaps(gapModel, modelDB, ''algorithm'', ''fastLP'', ''verbose'', false);');
            testCase.verifyClass(newModel, 'struct');
            testCase.verifyEqual(exitFlag, 1);
        end

        function fillGapsDispatchesGapfillMILP(testCase)
            testCase.assumeMILPSolver();
            modelDB = testCase.model; modelDB.id = 'DB';
            gapModel = removeReactions(modelDB, modelDB.rxns(1:5));
            gapModel.id = 'gapModel';
            evalc('[~,~,addedRxns,newModel,exitFlag] = fillGaps(gapModel, modelDB, ''algorithm'', ''gapfillMILP'', ''verbose'', false);');
            testCase.verifyClass(newModel, 'struct');
            testCase.verifyEqual(exitFlag, 1);
        end

    end

    methods (Access = private)
        function m = chainModel(~)
            % A -> B -> C, no exchanges; the caller adds the ones it needs.
            m = struct();
            m.rxns = {'R1';'R2'}; m.rxnNames = {'R1';'R2'};
            m.mets = {'A';'B';'C'}; m.metNames = {'A';'B';'C'}; m.metComps = [1;1;1];
            m.comps = {'c'}; m.compNames = {'c'};
            m.S = sparse([-1 0; 1 -1; 0 1]);
            m.lb = [0;0]; m.ub = [1000;1000]; m.rev = [0;0]; m.c = [0;1];
            m.b = [0;0;0];
            m.genes = {}; m.grRules = {'';''}; m.rxnGeneMat = sparse(2,0);
        end
    end
end
