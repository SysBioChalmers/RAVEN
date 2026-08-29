classdef tSyntax < RavenTestCase
% tSyntax  Every RAVEN source file must parse.
%
%   A MATLAB file with unbalanced block keywords still sits in the repository
%   looking perfectly ordinary: the error only surfaces when something calls
%   it, so a file behind an unconditional skip (e.g. assumeFail) can carry a
%   broken parse indefinitely without any test catching it.
%
%   This test is the cheap backstop: it needs no data, no solver and no
%   network, and it runs over the entire source tree in seconds.

    methods (Test)
        function allSourceFilesParse(testCase)
            files = dir(fullfile(testCase.ravenRoot,'**','*.m'));
            paths = fullfile({files.folder}, {files.name});

            % software/ is vendored third-party code (GLPKmex, libSBML). It is
            % not ours to fix, so a syntax complaint there is not a RAVEN bug.
            paths = paths(~contains(paths, [filesep 'software' filesep]));

            testCase.assertNotEmpty(paths, 'Found no .m files to check.');

            msgs = checkcode(paths{:}, '-id');
            if ~iscell(msgs)
                msgs = {msgs};   % checkcode returns a bare struct for one file
            end

            offenders = {};
            for i = 1:numel(paths)
                if isempty(msgs{i})
                    continue
                end
                % SYNER is checkcode's parse-error identifier. Only syntax is
                % checked here; style warnings are deliberately ignored so the
                % test stays a hard gate rather than a lint backlog.
                bad = strcmp({msgs{i}.id}, 'SYNER');
                for k = find(bad)
                    offenders{end+1} = sprintf('%s (line %d): %s', ...
                        strrep(paths{i}, testCase.ravenRoot, ''), ...
                        msgs{i}(k).line, msgs{i}(k).message); %#ok<AGROW>
                end
            end

            testCase.verifyEmpty(offenders, ...
                sprintf('Files with parse errors:\n%s', strjoin(offenders, newline)));
        end
    end
end
