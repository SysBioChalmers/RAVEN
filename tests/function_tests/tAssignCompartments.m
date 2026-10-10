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
    end
end
