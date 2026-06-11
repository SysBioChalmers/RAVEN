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
