classdef tBinaries < matlab.unittest.TestCase
% tBinaries  On-demand download of RAVEN's external binaries from raven-data.
%
%   Validates downloadRavenBinaries on the current platform: each tool is
%   fetched into software/<tool>/, the executable is present, and (for the CLI
%   tools) it actually runs. On macOS, running the downloaded binary also
%   confirms the Gatekeeper quarantine attribute was cleared.
%
%   Exercised on Linux and macOS by the "binary-tests" CI job. Requires network
%   access to https://github.com/SysBioChalmers/raven-data .

    methods (Test)
        function blast(testCase)
            testCase.fetchAndRun('blast+', 'blastp', '-version');
            testCase.verifyTrue(isfile(testCase.binPath('blast+','makeblastdb')), ...
                'makeblastdb not found after download');
        end

        function diamond(testCase)
            testCase.fetchAndRun('diamond', 'diamond', 'version');
        end

        function hmmer(testCase)
            testCase.fetchAndRun('hmmer', 'hmmsearch', '-h');
        end
    end

    methods (Access = private)
        function p = binPath(~, tool, name)
            if ispc; ext = '.exe'; elseif ismac; ext = '.mac'; else; ext = ''; end
            p = fullfile(findRAVENroot(),'software',tool,[name ext]);
        end

        function fetchAndRun(testCase, tool, name, versionArg)
            % Download the tool, then confirm the executable exists and runs.
            downloadRavenBinaries({tool});
            exe = testCase.binPath(tool, name);
            testCase.verifyTrue(isfile(exe), ...
                sprintf('%s not found after download', name));
            [status, ~] = system(['"' exe '" ' versionArg]);
            testCase.verifyEqual(status, 0, ...
                sprintf('%s did not execute (exit status %d)', name, status));
        end
    end
end
