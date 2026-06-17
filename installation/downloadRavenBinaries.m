function downloadRavenBinaries(tools)
% downloadRavenBinaries
%   Download RAVEN's external command-line binaries on demand from the
%   raven-data repository into RAVEN's software/ directory. These binaries are
%   no longer committed to the RAVEN git repository; only the binary matching
%   the current platform is fetched, and only when needed.
%
%   Input:
%   tools   cell array, subset of {'blast+','diamond','hmmer','WoLFPSORT'} to
%           download. Defaults to all of them. Tools whose target executable is
%           already present are skipped.
%
%   This requires an internet connection. For an offline / air-gapped install,
%   download the "*-binaries" RAVEN release (which bundles all binaries) and
%   unzip it into the RAVEN software/ directory instead.
%
%   The binaries are hosted as release assets in
%   https://github.com/SysBioChalmers/raven-data.
%
%   Usage: downloadRavenBinaries({'blast+','diamond'})
%
% NOTE: lazily invoked by getBlast/getDiamond/getKEGGModelForOrganism/
% getWoLFScores and offered by checkInstallation when a binary is missing.

if nargin < 1 || isempty(tools)
    tools = {'blast+','diamond','hmmer','WoLFPSORT'};
elseif ischar(tools)
    tools = {tools};
end

ravenDir = findRAVENroot();
base     = 'https://github.com/SysBioChalmers/raven-data/releases/download';

% raven-data platform key for the per-platform binary ZIPs.
if ispc
    plat = 'windows-x86_64';
elseif ismac
    plat = 'macos-arm64';
else
    plat = 'linux-x86_64';
end

for i = 1:numel(tools)
    tool = tools{i};
    switch tool
        case 'blast+'
            tag = 'blast-2.17.0';
            asset = ['blast-2.17.0-' plat '.zip'];
            destDir = fullfile(ravenDir,'software','blast+');
            sentinel = fullfile(destDir,'blastp');           % bare name on every OS (.exe on Windows)
            execs = {'blastp','makeblastdb'};
        case 'diamond'
            tag = 'diamond-2.1.17';
            asset = ['diamond-2.1.17-' plat '.zip'];
            destDir = fullfile(ravenDir,'software','diamond');
            sentinel = fullfile(destDir,'diamond');
            execs = {'diamond'};
        case 'hmmer'
            if ispc
                tag = 'hmmer-3.3.2';                          % native Windows build (Cygwin); reads HMMER3/f
                asset = ['hmmer-3.3.2-' plat '.zip'];
            else
                tag = 'hmmer-3.4.0';
                asset = ['hmmer-3.4.0-' plat '.zip'];
            end
            destDir = fullfile(ravenDir,'software','hmmer');
            sentinel = fullfile(destDir,'hmmsearch');
            execs = {'hmmsearch'};
        case 'WoLFPSORT'
            tag = 'wolfpsort-0.2';
            asset = 'wolfpsort-0.2.zip';                     % zip contains the WoLFPSORT/ tree
            destDir = fullfile(ravenDir,'software');
            sentinel = fullfile(destDir,'WoLFPSORT','bin','wolfPredict');
            execs = {};                                      % whole bin/ tree (chmod'd below)
        otherwise
            error('downloadRavenBinaries: unknown tool ''%s''.',tool);
    end

    % macOS uses a .mac suffix and Windows a .exe suffix (Linux: bare name).
    if ~strcmp(tool,'WoLFPSORT')
        if ispc
            sentinel = [sentinel '.exe']; %#ok<AGROW>
        elseif ismac
            sentinel = [sentinel '.mac']; %#ok<AGROW>
        end
    end
    if exist(sentinel,'file')
        continue;   % already provisioned
    end

    url = [base '/' tag '/' asset];
    zipPath = fullfile(ravenDir,'software',asset);
    fprintf('Downloading %s from raven-data ...\n',tool);
    try
        websave(zipPath,url);
    catch
        error(['Failed to download %s from %s\n' ...
               'Check your internet connection, or fetch the offline ' ...
               '"*-binaries" RAVEN release.'],tool,url);
    end
    unzip(zipPath,destDir);
    delete(zipPath);

    % macOS uses the .mac suffix throughout RAVEN, but the raven-data macOS ZIP
    % ships bare names, so rename the extracted executables to match.
    if ismac && ~strcmp(tool,'WoLFPSORT')
        for k = 1:numel(execs)
            bare = fullfile(destDir,execs{k});
            if isfile(bare)
                movefile(bare, fullfile(destDir,[execs{k} '.mac']));
            end
        end
    end

    % MATLAB's unzip does not preserve Unix execute permissions; restore them.
    if isunix
        if strcmp(tool,'WoLFPSORT')
            attrTarget = fullfile(destDir,'WoLFPSORT');
            system(['chmod -R +x "' fullfile(attrTarget,'bin') '"']);
        else
            attrTarget = destDir;
            if ismac; sfx = '.mac'; else; sfx = ''; end
            for k = 1:numel(execs)
                system(['chmod +x "' fullfile(destDir,[execs{k} sfx]) '"']);
            end
        end
        % Downloaded binaries are quarantined by macOS Gatekeeper and will not
        % execute until the com.apple.quarantine attribute is cleared.
        if ismac
            system(['xattr -dr com.apple.quarantine "' attrTarget '"']);
        end
    end
end
end
