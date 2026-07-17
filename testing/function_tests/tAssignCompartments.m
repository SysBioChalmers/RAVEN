classdef tAssignCompartments < RavenTestCase
% tAssignCompartments  Tests for the functionality-constrained compartment-assignment MILP.
%
%   Uses a tiny single-compartment toy: EX_A (uptake) -> r1: A->B -> bio: B-> (objective), with
%   gene g1 (on r1) scored higher for a hypothetical compartment 'm' than for the cytosol 'c'.
%   The decisive checks are that the functionality constraint can override the localization score
%   and that the returned model stays able to grow. Requires a MILP solver.

    methods (Test)

        function functionalityOverridesScore(testCase)
            testCase.assumeMILPSolver();
            model = tAssignCompartments.toy();
            GSS = tAssignCompartments.gss();
            % B is non-transportable, so placing r1 in 'm' would strand B where biomass (in 'c')
            % cannot reach it -> r1 must stay in 'c', against g1's higher 'm' score.
            evalc(['[oModel, place, ~, eFlag] = assignCompartments(model, GSS, {''r1''}, ' ...
                   '''defaultCompartment'', ''c'', ''transportable'', {}, ''verbose'', false);']);
            testCase.verifyEqual(eFlag, 1);
            testCase.verifyEqual(place.compartment{1}, 'c');
            sol = solveLP(oModel);
            testCase.verifyGreaterThan(sol.f, 1e-6);          % applied model still grows
        end

        function scoreRespectedWhenFunctional(testCase)
            testCase.assumeMILPSolver();
            model = tAssignCompartments.toy();
            GSS = tAssignCompartments.gss();
            % Transports allowed and cheap -> r1 can follow its score to 'm' and stay functional.
            evalc(['[oModel, place, trans, eFlag] = assignCompartments(model, GSS, {''r1''}, ' ...
                   '''defaultCompartment'', ''c'', ''transportCost'', 0.1, ''verbose'', false);']);
            testCase.verifyEqual(eFlag, 1);
            testCase.verifyEqual(place.compartment{1}, 'm');
            testCase.verifyGreaterThan(numel(trans.mets), 0); % transports were added
            sol = solveLP(oModel);
            testCase.verifyGreaterThan(sol.f, 1e-6);
        end

        function infeasibleGrowthFloorReported(testCase)
            testCase.assumeMILPSolver();
            model = tAssignCompartments.toy();
            GSS = tAssignCompartments.gss();
            evalc(['[~, ~, ~, eFlag] = assignCompartments(model, GSS, {''r1''}, ' ...
                   '''defaultCompartment'', ''c'', ''minGrowth'', 1e6, ''verbose'', false);']);
            testCase.verifyEqual(eFlag, -1);                  % growth floor unreachable
        end

        function certificationReportsRealGrowth(testCase)
            % A certified placement reports certified=true and the achieved
            % growth; an unreachable floor reports the real (short) growth
            % rather than hiding it behind an optimal-MILP status.
            testCase.assumeMILPSolver();
            model = tAssignCompartments.toy();
            GSS = tAssignCompartments.gss();
            evalc(['[~, ~, ~, ok, rep] = assignCompartments(model, GSS, {''r1''}, ' ...
                   '''defaultCompartment'', ''c'', ''transportable'', {}, ''verbose'', false);']);
            testCase.verifyEqual(ok, 1);
            testCase.verifyTrue(rep.certified);
            testCase.verifyEqual(rep.status, 'certified');
            testCase.verifyGreaterThan(rep.growths.primary, 1e-6);

            evalc(['[~, ~, ~, bad, badRep] = assignCompartments(model, GSS, {''r1''}, ' ...
                   '''defaultCompartment'', ''c'', ''minGrowth'', 1e6, ''verbose'', false);']);
            testCase.verifyEqual(bad, -1);
            testCase.verifyFalse(badRep.certified);
            testCase.verifyEqual(badRep.status, 'uncertified');
            % the real growth is small and reported, not concealed
            testCase.verifyGreaterThan(badRep.growths.primary, 1e-6);
            testCase.verifyLessThan(badRep.growths.primary, 1e6);
        end

        function confinementColocatesMovableReactions(testCase)
            % Two movable reactions share a non-transportable intermediate X.
            % Their scores pull them to different compartments, but X cannot
            % be transported, so the confinement repair must co-locate them.
            testCase.assumeMILPSolver();
            model = tAssignCompartments.chainToy();
            GSS = tAssignCompartments.chainGss();
            % X non-transportable (S/P transportable, keyed by metabolite
            % name), so the biomass path is not what forces co-location.
            evalc(['[~, place, ~, ok] = assignCompartments(model, GSS, {''r1'';''r2''}, ' ...
                   '''defaultCompartment'', ''c'', ''transportable'', {''S'';''P''}, ''verbose'', false);']);
            testCase.verifyEqual(ok, 1);
            c1 = place.compartment{strcmp(place.rxns,'r1')};
            c2 = place.compartment{strcmp(place.rxns,'r2')};
            testCase.verifyEqual(c1, c2);                     % co-located

            % With X transportable too, nothing forces co-location and the
            % dominant scores split the two reactions across compartments.
            evalc(['[~, sp] = assignCompartments(model, GSS, {''r1'';''r2''}, ' ...
                   '''defaultCompartment'', ''c'', ''transportable'', {''S'';''X'';''P''}, ''verbose'', false);']);
            testCase.verifyNotEqual(sp.compartment{strcmp(sp.rxns,'r1')}, ...
                                    sp.compartment{strcmp(sp.rxns,'r2')});
        end

    end

    methods (Static)
        function model = toy()
            model.id='toy'; model.comps={'c'}; model.compNames={'cytoplasm'};
            model.mets={'A_c';'B_c'}; model.metNames={'A';'B'}; model.metComps=[1;1];
            model.S=sparse([-1 -1 0; 0 1 -1]);     % [A;B] x [EX_A r1 bio]
            model.rxns={'EX_A';'r1';'bio'}; model.rxnNames=model.rxns;
            model.lb=[-10;0;0]; model.ub=[1000;1000;1000]; model.rev=[1;0;0];
            model.c=[0;0;1]; model.b=zeros(2,1);
            model.genes={'g1'}; model.grRules={'';'g1';''}; model.rxnGeneMat=sparse([0;1;0]);
        end
        function GSS = gss()
            GSS.genes={'g1'}; GSS.compartments={'c';'m'}; GSS.scores=[0.4 0.9];
        end

        function model = chainToy()
            % EX_S -> S -> [r1] -> X -> [r2] -> P -> bio.  X is the shared,
            % non-transportable intermediate of the two movable reactions.
            model.id='chain'; model.comps={'c'}; model.compNames={'cytoplasm'};
            model.mets={'S_c';'X_c';'P_c'}; model.metNames={'S';'X';'P'}; model.metComps=[1;1;1];
            %              EX_S  r1    r2    bio
            model.S=sparse([ 1   -1    0     0;    % S
                             0    1   -1     0;    % X
                             0    0    1    -1]);  % P
            model.rxns={'EX_S';'r1';'r2';'bio'}; model.rxnNames=model.rxns;
            model.lb=[-10;0;0;0]; model.ub=[1000;1000;1000;1000]; model.rev=[1;0;0;0];
            model.c=[0;0;0;1]; model.b=zeros(3,1);
            model.genes={'g1';'g2'}; model.grRules={'';'g1';'g2';''};
            model.rxnGeneMat=sparse([0 0; 1 0; 0 1; 0 0]);
        end
        function GSS = chainGss()
            % Each gene has a single dominant compartment (the other score is
            % below the 0.5 multi-compartment penalty), so placement is
            % determined without a tie-break: g1 -> 'c', g2 -> 'm'. When X is
            % transportable they split; when X is non-transportable the shared
            % pool co-locates them, in 'c' (combined 0.9+0.3 > 0.1+0.9).
            GSS.genes={'g1';'g2'}; GSS.compartments={'c';'m'};
            GSS.scores=[0.9 0.1; 0.3 0.9];
        end
    end
end
