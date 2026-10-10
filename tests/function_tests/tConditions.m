classdef tConditions < RavenTestCase
% tConditions  Tests for the condition-application functions in conditions/.

    methods (Test)

        function applyConditionSetsBounds(testCase)
            rxn = testCase.model.rxns{1};
            condition.bounds = {struct('rxn', rxn, 'lb', 0, 'ub', 0)};
            evalc('m2 = applyCondition(testCase.model, condition);');
            idx = strcmp(m2.rxns, rxn);
            testCase.verifyEqual(m2.ub(idx), 0);
            testCase.verifyEqual(m2.lb(idx), 0);
        end

        function applyConditionResetExchangesZeroesUptake(testCase)
            % The prelude reopens exchanges to (0, 1000), i.e. no uptake.
            m = testCase.model;
            [~, exchangeRxns] = getExchangeRxns(m);
            m.lb(exchangeRxns) = -37;
            m.ub(exchangeRxns) = 42;
            condition.prelude.reset_exchanges = 'both';
            evalc('m2 = applyCondition(m, condition);');
            testCase.verifyEqual(unique(m2.lb(exchangeRxns)), 0);
            testCase.verifyEqual(unique(m2.ub(exchangeRxns)), 1000);
        end

        function applyConditionRemoveMetZeroesCoefficient(testCase)
            % remove_mets drops a metabolite from the cofactor pseudoreaction.
            m = tConditions.cofactorModel();
            condition.cofactor_pseudoreaction = struct('rxn_id', 'R_cofactor', ...
                'remove_mets', {{struct('met', 'nadh_c')}});
            evalc('m2 = applyCondition(m, condition);');
            testCase.verifyEqual(full(m2.S(strcmp(m2.mets,'nadh_c'), 1)), 0);
        end

        function applyConditionRecomputesChargeBalance(testCase)
            % charge_balance_met takes whatever coefficient makes the
            % pseudoreaction charge neutral, recomputed after the removals.
            m = tConditions.cofactorModel();
            condition.cofactor_pseudoreaction = struct('rxn_id', 'R_cofactor', ...
                'remove_mets', {{struct('met', 'nadh_c')}}, ...
                'charge_balance_met', 'h_c');
            evalc('m2 = applyCondition(m, condition);');
            % Left after removing NADH: atp_c (-1, charge -4) and amp_c
            % (+1, charge -2), a residual of +2, so H+ must enter at -2.
            testCase.verifyEqual(full(m2.S(strcmp(m2.mets,'h_c'), 1)), -2);
            % which is to say the reaction ends up charge neutral
            idx = find(m2.S(:,1));
            testCase.verifyEqual(sum(full(m2.S(idx,1)) .* m2.metCharges(idx)), 0);
        end

        function applyConditionBiomassDeltaAddsToCoefficient(testCase)
            % A stoichiometry delta adds to the existing coefficient rather
            % than replacing it.
            m = tConditions.cofactorModel();
            m.rxns{2} = 'R_biomass'; m.rxnNames{2} = 'R_biomass';
            before = full(m.S(strcmp(m.mets,'atp_c'), 2));
            condition.biomass_stoichiometry_delta = struct('rxn_id', 'R_biomass', ...
                'add', {{struct('met', 'atp_c', 'coef', -0.5)}});
            evalc('m2 = applyCondition(m, condition);');
            testCase.verifyEqual(full(m2.S(strcmp(m2.mets,'atp_c'), 2)), before - 0.5);
        end

    end

    methods (Static, Access = private)

        function m = cofactorModel()
            % R_cofactor consumes ATP and NADH; R_other consumes ATP.
            m = struct();
            m.id        = 'cofac';
            m.rxns      = {'R_cofactor';'R_other'};
            m.rxnNames  = {'R_cofactor';'R_other'};
            m.mets      = {'atp_c';'nadh_c';'amp_c';'h_c'};
            m.metNames  = {'ATP';'NADH';'AMP';'H+'};
            m.metComps  = [1;1;1;1];
            m.metCharges= [-4;-2;-2;1];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            %              R_cofactor R_other
            m.S = sparse([  -1        -1;    % atp_c
                            -1         0;    % nadh_c
                             1         0;    % amp_c
                             0         0]);  % h_c
            m.lb = [0;0]; m.ub = [1000;1000]; m.rev = [0;0]; m.c = [0;0];
            m.b = zeros(4,1);
            m.genes = {}; m.grRules = {'';''}; m.rxnGeneMat = sparse(2,0);
            m.metFormulas = {'C10';'C21';'C10';'H'};
        end

    end
end
