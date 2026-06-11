classdef tLocalization < RavenTestCase
% tLocalization  Tests for the subcellular-localization functions in localization/.
%
%   This module is driven by the external WoLF PSORT predictor (Linux-only)
%   and by predictor output files / gene-score structures that are produced
%   by it. Without that pipeline the functions cannot run, so each is guarded
%   and skips in this environment.

    methods (Test)

        function getWoLFScoresNeedsPredictor(testCase)
            testCase.assumeFail('Requires the external WoLF PSORT predictor (Linux).');
        end

        function parseScoresNeedsPredictorOutput(testCase)
            testCase.assumeFail('Requires a predictor output file.');
        end

        function mapCompartmentsNeedsGeneScoreStructure(testCase)
            testCase.assumeFail('Requires a gene-score structure from getWoLFScores/parseScores.');
        end

        function predictLocalizationNeedsGeneScoreStructure(testCase)
            testCase.assumeFail('Requires a gene-score structure from getWoLFScores/parseScores.');
        end

    end
end
