classdef tInstallation < RavenTestCase
% tInstallation  Tests for the setup/installation helpers in installation/.
%
%   addRavenToUserPath, removeRavenFromPath and updateDocumentation have
%   destructive / irreversible side effects (rewrite the user's startup.m,
%   remove RAVEN from the live path, regenerate the whole doc tree). They are
%   therefore only checked for existence here, not executed.

    methods (Test)

        function findRAVENrootReturnsInstallFolder(testCase)
            p = findRAVENroot();
            testCase.verifyTrue(isfolder(p));
            testCase.verifyTrue(exist(fullfile(p,'installation','findRAVENroot.m'),'file')==2);
        end

        function checkInstallationReturnsVersion(testCase)
            [~, currVer] = evalc('checkInstallation(false, false)');
            testCase.verifyNotEmpty(currVer);
        end

        function checkFunctionUniquenessRuns(testCase)
            [~, status] = evalc('checkFunctionUniqueness()');
            testCase.verifyTrue(islogical(status) || isnumeric(status));
        end

        function addRavenToUserPathExists(testCase)
            testCase.verifyEqual(exist('addRavenToUserPath','file'), 2);
        end

        function removeRavenFromPathExists(testCase)
            testCase.verifyEqual(exist('removeRavenFromPath','file'), 2);
        end

        function updateDocumentationExists(testCase)
            testCase.verifyEqual(exist('updateDocumentation','file'), 2);
        end

    end
end
