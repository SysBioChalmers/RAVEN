classdef tBiomass < RavenTestCase
% tBiomass  Tests for the biomass-composition functions in biomass/.
%
%   The textbook E. coli model has a single lumped biomass reaction rather
%   than per-component pseudoreactions, so component-scaling functions are
%   tested on their missing-pseudoreaction error path, and fitParameters
%   (which needs quadprog) is guarded on the Optimization Toolbox.

    methods (Test)

        function getBiomassFractionsHandlesMissingComponents(testCase)
            cfg = testCase.biomassConfig();
            fractions = getBiomassFractions(testCase.model, cfg);
            testCase.verifyTrue(isfield(fractions, 'total'));
        end

        function setGAMReturnsModel(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            evalc(['m2 = setGAM(testCase.model, 50, biomassRxn, ' ...
                '{''ATP'',''ADP'',''H2O'',''H+'',''Phosphate''});']);
            testCase.verifyClass(m2, 'struct');
        end

        function scaleBiomassFractionErrorsOnMissingComponent(testCase)
            cfg = testCase.biomassConfig();
            testCase.verifyError(@() scaleBiomassFraction(testCase.model, cfg, ...
                'protein', 0.5), ?MException);
        end

        function scaleBiomassPseudoreactionErrorsOnMissingComponent(testCase)
            cfg = testCase.biomassConfig();
            testCase.verifyError(@() scaleBiomassPseudoreaction(testCase.model, cfg, ...
                'protein', 0.9), ?MException);
        end

        function scaleBiomassPseudoreactionRebalancesProton(testCase)
            % Rescaling the substrates rebalances H+ so the pseudoreaction
            % stays charge neutral.
            [m, cfg] = tBiomass.protonToyModel();
            out = scaleBiomassPseudoreaction(m, cfg, 'protein', 0.5);
            % x_c: -0.5 * charge -2 = +1, so H+ must come in at -1
            testCase.verifyEqual(full(out.S(strcmp(out.mets,'x_c'), 1)), -0.5, 'AbsTol', 1e-9);
            testCase.verifyEqual(full(out.S(strcmp(out.mets,'h_c'), 1)), -1, 'AbsTol', 1e-9);
        end

        function scaleBiomassPseudoreactionRefusesUnknownCharge(testCase)
            % An unset charge on a participating metabolite leaves the charge
            % balance unknown. Treating it as neutral would silently write the
            % wrong proton coefficient, so the rescale must refuse instead.
            [m, cfg] = tBiomass.protonToyModel();
            m.metCharges(strcmp(m.mets,'x_c')) = NaN;
            testCase.verifyError(@() scaleBiomassPseudoreaction(m, cfg, 'protein', 0.5), ...
                'scaleBiomassPseudoreaction:unknownCharge');
        end

        function fitParametersRunsWhenQuadprogAvailable(testCase)
            testCase.assumeDependency(exist('quadprog','file')==2, ...
                'Optimization Toolbox (quadprog)');
            % Minimal single-parameter fit of an exchange flux.
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            pos.position = 1; pos.value = 0;
            evalc(['p = fitParameters(testCase.model, {biomassRxn}, 0.5, ' ...
                '{biomassRxn}, 0.5, pos);']);
            testCase.verifyNotEmpty(p);
        end

    end

    methods (Static, Access = private)
        function [m, cfg] = protonToyModel()
            % One pseudoreaction: x_c -> protein, with h_c available to carry
            % the charge balance.
            m = struct();
            m.id        = 'toy';
            m.rxns      = {'R1'};
            m.rxnNames  = {'protein pseudoreaction'};
            m.mets      = {'x_c';'h_c';'prot_c'};
            m.metNames  = {'x';'h';'protein'};
            m.metComps  = [1;1;1];
            m.metCharges= [-2;1;0];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-1;0;1]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 1; m.b = zeros(3,1);
            m.grRules = {''}; m.genes = {}; m.rxnGeneMat = sparse(1,0);

            cfg = struct();
            cfg.biomass_rxn = 'R1';
            cfg.proton_met  = 'h_c';
            cfg.components  = {struct('name','protein', ...
                'pseudoreaction_name','protein pseudoreaction', ...
                'mass_strategy','mw')};
        end
    end

    methods (Access = private)
        function cfg = biomassConfig(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            cfg.biomass_rxn = biomassRxn;
            cfg.proton_met = 'h_c';
            cfg.components = {struct('name','protein', ...
                'pseudoreaction_name','no_such_pseudoreaction','mass_strategy','mw')};
        end
    end
end
