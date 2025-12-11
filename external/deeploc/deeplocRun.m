function [csvFilePath, pollURL] = deeplocRun_old(fastaFilePath, varargin)
% deeplocRun
%   Submits protein sequences from a FASTA file to DeepLoc webserver
%   and downloads the resulting CSV file with localization predictions.
%
% Input:
%   fastaFilePath     path to FASTA file containing protein sequences
%   'Version'         DeepLoc version: '2.1' (default), '2.0', or '1.0'
%                     (optional)
%   'Mode'            prediction mode: 'high-quality' (default) or 'fast'
%                     (optional, for DeepLoc 2.0/2.1 only)
%   'Figures'         true if figures should be requested from DeepLoc
%                     (optional, default false, for DeepLoc 2.0/2.1 only)
%   'Encoding'        protein encoding: 'profiles' (default) or 'blosum62'
%                     (optional, for DeepLoc 1.0 only)
%   'Email'           email address for notification when job completes
%                     (optional, default: no email)
%   'OutputDir'       directory to save CSV file (optional, default: same
%                     directory as FASTA file)
%   'Timeout'         HTTP request timeout in seconds (optional, default 900)
%   'Verbose'         true if progress messages should be printed (optional,
%                     default true)
%
% Output:
%   csvFilePath       full path to downloaded CSV file with DeepLoc results
%   pollURL           URL to check job status on DeepLoc webserver
%
% Usage:
%   [csvFile, pollURL] = deeplocRun('proteins.fasta');
%   [csvFile, pollURL] = deeplocRun('proteins.fasta', 'Version', '2.0');
%   [csvFile, pollURL] = deeplocRun('proteins.fasta', 'Version', '1.0', 'Encoding', 'blosum62');
%   [csvFile, pollURL] = deeplocRun('proteins.fasta', 'Mode', 'fast', ...
%       'Figures', true, 'Email', 'user@example.com');
%
% NOTE: This function requires internet connection to access DeepLoc-2.1
% webserver. The CSV file format must be compatible with parseScores.m
% using the 'deeploc' option.
%
% The workflow is:
%   1. deeplocRun() → CSV file
%   2. parseScores(CSV, 'deeploc') → GSS structure
%   3. predictLocalization(model, GSS, ...) → compartmentalized model

% Valid DeepLoc versions
validVersions = {'2.1', '2.0', '1.0'};

% Parse optional parameters
p = inputParser;
addParameter(p, 'Version', '2.1', @(x) ischar(x) || isstring(x));
addParameter(p, 'Mode', 'high-quality', @(x) ischar(x) || isstring(x));
addParameter(p, 'Figures', false, @islogical);
addParameter(p, 'Encoding', 'profiles', @(x) ischar(x) || isstring(x));
addParameter(p, 'Email', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x) || isempty(x));  % Will be set to FASTA directory if empty
addParameter(p, 'Timeout', 900, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});

% Get version from parser (allows override via named parameter)
version = char(p.Results.Version);

% Validate version
if ~ismember(version, validVersions)
    EM = sprintf('Version ''%s'' is not supported. Valid versions are: %s', version, strjoin(validVersions, ', '));
    dispEM(EM, true);
end

% Construct baseURL from version
baseURL = sprintf('https://services.healthtech.dtu.dk/services/DeepLoc-%s', version);

mode = char(p.Results.Mode);
figures = p.Results.Figures;
encoding = char(p.Results.Encoding);
email = char(p.Results.Email);
outputDir = char(p.Results.OutputDir);
timeout = p.Results.Timeout;
verbose = p.Results.Verbose;

% Normalize inputs
fastaFilePath = convertCharArray(fastaFilePath);
if numel(fastaFilePath) > 1
    EM = 'Only one FASTA file can be submitted at a time';
    dispEM(EM, true);
end
fastaFilePath = fastaFilePath{1};

% Validate FASTA file exists
fastaFilePath = checkFileExistence(fastaFilePath, 1);

% Set default output directory to FASTA file directory if not specified
if isempty(outputDir)
    [fastaDir, ~, ~] = fileparts(fastaFilePath);
    outputDir = fastaDir;
end

% Validate output directory
if ~isfolder(outputDir)
    try
        mkdir(outputDir);
    catch
        EM = sprintf('Cannot create output directory: %s', outputDir);
        dispEM(EM, true);
    end
end

% Read and validate FASTA file
if verbose
    fprintf('Reading FASTA file: %s\n', fastaFilePath);
end

% Try to read FASTA file to validate format and check sequences
try
    % Try using fastaread if Bioinformatics Toolbox is available
    if exist('fastaread', 'file') == 2
        sequences = fastaread(fastaFilePath);
        if isempty(sequences)
            EM = 'FASTA file appears to be empty';
            dispEM(EM, true);
        end
        numSequences = numel(sequences);
        if verbose
            fprintf('  Found %d sequence(s)\n', numSequences);
        end
        
        % Check sequence lengths and count
        shortSequences = false;
        longSequences = false;
        for i = 1:numSequences
            seqLength = length(sequences(i).Sequence);
            if seqLength < 10
                shortSequences = true;
            end
            if seqLength > 6000
                longSequences = true;
            end
        end
        
        % Issue warnings for short/long sequences
        if shortSequences || longSequences
            warningMsg = 'Warning: ';
            if shortSequences && longSequences
                warningMsg = [warningMsg 'Amino acid sequences shorter than 10 amino acids and longer than 6000 amino acids have been detected. '];
            elseif shortSequences
                warningMsg = [warningMsg 'Amino acid sequences shorter than 10 amino acids have been detected. '];
            else
                warningMsg = [warningMsg 'Amino acid sequences longer than 6000 amino acids have been detected. '];
            end
            warningMsg = [warningMsg 'DeepLoc predictions may be inaccurate for those sequences.'];
            if verbose
                fprintf('  %s\n', warningMsg);
            else
                dispEM(warningMsg, false);
            end
        end
    else
        % Simple validation: check if file contains '>' characters
        fid = fopen(fastaFilePath, 'r');
        if fid == -1
            EM = sprintf('Cannot open FASTA file: %s', fastaFilePath);
            dispEM(EM, true);
        end
        content = fread(fid, '*char')';
        fclose(fid);
        if ~contains(content, '>')
            EM = 'FASTA file does not appear to contain valid FASTA format (no headers starting with ''>'')';
            dispEM(EM, false);
        end
        
        % Count sequences by counting '>' characters
        numSequences = sum(content == '>');
        if verbose
            fprintf('  Found %d sequence(s)\n', numSequences);
        end
        
        % For basic validation without fastaread, we can't check sequence lengths easily
        % So we'll skip the length warnings in this case
    end
catch ME
    EM = sprintf('Error reading FASTA file: %s', ME.message);
    dispEM(EM, true);
end

% Extract version from baseURL (e.g., "2.1" from "https://.../DeepLoc-2.1")
% This must be done early to determine which parameters are valid
versionMatch = regexp(baseURL, 'DeepLoc-([0-9]+\.[0-9]+)', 'tokens', 'once');
if isempty(versionMatch)
    EM = sprintf('Could not extract version from baseURL: %s. Expected format: .../DeepLoc-X.Y', baseURL);
    dispEM(EM, true);
end
deeplocVersion = versionMatch{1};
versionMajor = str2double(regexp(deeplocVersion, '^([0-9]+)', 'tokens', 'once'));

% Validate parameters based on version
if versionMajor >= 2
    % DeepLoc 2.0/2.1: validate Mode
    if ~ismember(lower(mode), {'fast', 'high-quality'})
        EM = 'Mode must be either ''fast'' or ''high-quality'' (for DeepLoc 2.0/2.1)';
        dispEM(EM, true);
    end
else
    % DeepLoc 1.0: validate Encoding
    if ~ismember(lower(encoding), {'profiles', 'blosum62'})
        EM = 'Encoding must be either ''profiles'' or ''blosum62'' (for DeepLoc 1.0)';
        dispEM(EM, true);
    end
end

% Validate sequence count based on version and encoding
if exist('numSequences', 'var')
    maxSequences = 500;  % Default maximum
    if versionMajor < 2 && strcmpi(encoding, 'profiles')
        % DeepLoc 1.0 with Profiles encoding: maximum 50 sequences
        maxSequences = 50;
    end
    
    if numSequences > maxSequences
        if versionMajor < 2 && strcmpi(encoding, 'profiles')
            EM = sprintf('FASTA file contains %d sequences, but DeepLoc 1.0 with Profiles encoding allows a maximum of 50 sequences. Please reduce the number of sequences or use BLOSUM62 encoding (allows up to 500 sequences).', numSequences);
        else
            EM = sprintf('FASTA file contains %d sequences, but DeepLoc allows a maximum of %d sequences. Please reduce the number of sequences.', numSequences, maxSequences);
        end
        dispEM(EM, true);
    end
end

% Map parameters to DeepLoc values based on version
if versionMajor >= 2
    % DeepLoc 2.0/2.1: Map mode parameter to DeepLoc values
    if strcmpi(mode, 'fast')
        encodeMode = 'Fast';
    else
        encodeMode = 'Slow';  % 'high-quality' maps to 'Slow'
    end
    
    % Map figures parameter to DeepLoc format values
    if figures
        formatMode = 'long';
    else
        formatMode = 'short';
    end
else
    % DeepLoc 1.0: Map encoding parameter (values must be uppercase)
    if strcmpi(encoding, 'blosum62')
        encodingMode = 'BLOSUM62';
    else
        encodingMode = 'PROFILES';  % 'profiles' is default, uppercase for 1.0
    end
end

% Construct submission URL
urlParts = regexp(baseURL, '^(https?://[^/]+)', 'tokens', 'once');
if isempty(urlParts)
    EM = sprintf('Invalid baseURL format: %s', baseURL);
    dispEM(EM, true);
end
domainRoot = urlParts{1};
submitURL = [domainRoot '/cgi-bin/webface2.cgi'];

if verbose
    fprintf('Submitting to DeepLoc-%s...\n', deeplocVersion);
    if versionMajor >= 2
        fprintf('  Mode: %s (encode=%s)\n', mode, encodeMode);
        fprintf('  Format: %s\n', formatMode);
    else
        fprintf('  Encoding: %s\n', encodingMode);
    end
end

% ============================================================================
% STEP 1: Manually construct multipart/form-data body
% ============================================================================
if verbose
    fprintf('  Step 1: Preparing multipart form data...\n');
end

% Construct configfile path based on version (stable, from browser-captured request)
configfile = sprintf('/var/www/services/services/DeepLoc-%s/webface.cf', deeplocVersion);

% Read FASTA file as binary
fid = fopen(fastaFilePath, 'rb');
if fid == -1
    EM = sprintf('Cannot open FASTA file for reading: %s', fastaFilePath);
    dispEM(EM, true);
end
fastaFileData = fread(fid, '*uint8')';  % Ensure row vector
fclose(fid);
% Ensure fastaFileData is a row vector
if size(fastaFileData, 1) > 1
    fastaFileData = fastaFileData';
end

% Get filename for upload
[~, fastaFileName, fastaExt] = fileparts(fastaFilePath);
if isempty(fastaExt)
    fastaExt = '.faa';
end
fastaFileNameFull = [fastaFileName fastaExt];

% Generate random boundary string
boundary = ['----MATLABFormBoundary', char(randi([65 90], 1, 16))];

% Build multipart body manually
% Each part: --boundary\r\nContent-Disposition: ...\r\n\r\n<data>\r\n
bodyParts = {};

% Part 1: configfile
bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="configfile"\r\n');
bodyParts{end+1} = sprintf('\r\n');
bodyParts{end+1} = sprintf('%s\r\n', configfile);

% Part 2: Version-specific parameters (order matters for DeepLoc 1.0)
if versionMajor >= 2
    % DeepLoc 2.0/2.1: encode and format
    bodyParts{end+1} = sprintf('--%s\r\n', boundary);
    bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="encode"\r\n');
    bodyParts{end+1} = sprintf('\r\n');
    bodyParts{end+1} = sprintf('%s\r\n', encodeMode);
    
    bodyParts{end+1} = sprintf('--%s\r\n', boundary);
    bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="format"\r\n');
    bodyParts{end+1} = sprintf('\r\n');
    bodyParts{end+1} = sprintf('%s\r\n', formatMode);
else
    % DeepLoc 1.0: empty "fasta" field (for pasted sequences, must be before uploadfile)
    bodyParts{end+1} = sprintf('--%s\r\n', boundary);
    bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="fasta"\r\n');
    bodyParts{end+1} = sprintf('\r\n');
    bodyParts{end+1} = sprintf('\r\n');  % Empty field
end

% Part 3: uploadfile (binary file)
bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="uploadfile"; filename="%s"\r\n', fastaFileNameFull);
bodyParts{end+1} = sprintf('Content-Type: application/octet-stream\r\n');
bodyParts{end+1} = sprintf('\r\n');

% Closing boundary
bodyParts{end+1} = sprintf('--%s--\r\n', boundary);

% Convert all parts to uint8 and concatenate
bodyChar = [bodyParts{:}];
bodyUint8 = uint8(bodyChar);
% Ensure bodyUint8 is a row vector
if size(bodyUint8, 1) > 1
    bodyUint8 = bodyUint8';
end

% Insert binary file data before the closing boundary
% Find position of last boundary (before --boundary--)
closingBoundaryPos = length(bodyUint8) - length(sprintf('--%s--\r\n', boundary)) + 1;
% Ensure all parts are row vectors before concatenation
bodyUint8 = [bodyUint8(1:closingBoundaryPos-1), fastaFileData(:)', bodyUint8(closingBoundaryPos:end)];

% For DeepLoc 1.0, add encode field AFTER file data but BEFORE closing boundary
if versionMajor < 2
    % DeepLoc 1.0: encode field (same name as 2.0/2.1, but different values)
    % Values are "PROFILES" or "BLOSUM62" (uppercase)
    % Insert this right before the closing boundary
    encodeField = sprintf('--%s\r\nContent-Disposition: form-data; name="encode"\r\n\r\n%s\r\n', boundary, encodingMode);
    encodeFieldUint8 = uint8(encodeField);
    if size(encodeFieldUint8, 1) > 1
        encodeFieldUint8 = encodeFieldUint8';
    end
    % Insert encode field before closing boundary
    closingBoundaryPos = length(bodyUint8) - length(sprintf('--%s--\r\n', boundary)) + 1;
    bodyUint8 = [bodyUint8(1:closingBoundaryPos-1), encodeFieldUint8, bodyUint8(closingBoundaryPos:end)];
end

if verbose
    fprintf('    Data prepared for submission (%d bytes)\n', length(bodyUint8));
end

% ============================================================================
% STEP 2: Submit POST request and capture 302 redirect with jobid
% ============================================================================
if verbose
    fprintf('  Step 2: Submitting POST request to: %s\n', submitURL);
end

% Use matlab.net.http to send POST and capture 302 redirect headers
% (webwrite doesn't expose headers, so we need this for redirect handling)
import matlab.net.*
import matlab.net.http.*
import matlab.net.http.field.*
import matlab.net.http.io.*

% Create headers
% Use row vector concatenation to avoid dimension mismatch
contentTypeHeader = GenericField('Content-Type', sprintf('multipart/form-data; boundary=%s', boundary));
userAgentHeader = GenericField('User-Agent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36');
requestHeaders = [contentTypeHeader, userAgentHeader];  % Use comma (row vector) instead of semicolon

% Create request with body
% Convert uint8 to char (FASTA files are text, so this is safe)
% Wrap in MessageBody for proper handling
bodyChar = char(bodyUint8);
bodyMessage = MessageBody(bodyChar);
request = RequestMessage('POST', requestHeaders, bodyMessage);

% Set HTTP options to NOT follow redirects (we need to capture the 302)
httpOptions = HTTPOptions('MaxRedirects', 0, 'ConnectTimeout', timeout);

try
    % Send request
    [response, ~, ~] = request.send(URI(submitURL), httpOptions);
    
    if verbose
        fprintf('    POST request completed (status: %d)\n', response.StatusCode);
    end
    
    % Check for 302 redirect
    if response.StatusCode == 302
        % Extract Location header
        locationFields = response.getFields('Location');
        if isempty(locationFields)
            EM = 'Received 302 redirect but no Location header found';
            dispEM(EM, true);
        end
        
        % Get the Location header value (may be relative or absolute)
        locationURL = char(locationFields(1).Value);
        
        if verbose
            fprintf('    Location header value: %s\n', locationURL);
        end
        
        % Make absolute URL if relative
        if startsWith(locationURL, '/')
            locationURL = [domainRoot locationURL];
        elseif ~startsWith(locationURL, 'http')
            % If it doesn't start with http and isn't relative, prepend domain
            locationURL = [domainRoot '/' locationURL];
        end
        
        if verbose
            fprintf('    Full redirect URL: %s\n', locationURL);
        end
        
        % Parse jobid and wait value from Location URL
        % Format: /cgi-bin/webface2.cgi?jobid=XXXXXXXXXXXX&wait=20
        jobidPattern = '[?&]jobid=([^&]+)';
        jobidMatch = regexp(locationURL, jobidPattern, 'tokens', 'once');
        if isempty(jobidMatch)
            EM = sprintf('Could not extract jobid from redirect Location header. URL: %s', locationURL);
            dispEM(EM, true);
        end
        jobId = jobidMatch{1};
        
        waitPattern = '[?&]wait=(\d+)';
        waitMatch = regexp(locationURL, waitPattern, 'tokens', 'once');
        waitTime = 20;  % Default wait time
        if ~isempty(waitMatch)
            waitTime = str2double(waitMatch{1});
        end
        
        if verbose
            fprintf('    Extracted jobid: %s\n', jobId);
        end
        
    else
        % Not a redirect - check for errors
        % Extract response body for error checking
        if isa(response.Body, 'matlab.net.http.MessageBody')
            bodyData = response.Body.Data;
            if isa(bodyData, 'uint8')
                responseBody = char(bodyData');
            else
                responseBody = char(bodyData);
            end
        else
            responseBody = char(response.Body);
        end
        
        if contains(responseBody, 'Jobid not provided', 'IgnoreCase', true) || ...
           contains(responseBody, 'WebfaceConfigError', 'IgnoreCase', true)
            EM = sprintf('Server returned error (status %d): %s', response.StatusCode, responseBody);
            dispEM(EM, true);
        else
            EM = sprintf('Unexpected response status: %d (expected 302 redirect)', response.StatusCode);
            dispEM(EM, true);
        end
    end
    
catch ME
    EM = sprintf('Error submitting POST request to DeepLoc: %s', ME.message);
    dispEM(EM, true);
end

% Ensure jobId and waitTime variables exist
if ~exist('jobId', 'var') || isempty(jobId)
    EM = 'Could not extract jobid from server response';
    dispEM(EM, true);
end
if ~exist('waitTime', 'var')
    waitTime = 20;  % Default wait time if not extracted from redirect
end

% ============================================================================
% STEP 3: Poll for job completion using jobid
% ============================================================================
if verbose
    fprintf('  Step 3: Polling for job completion...\n');
end

% Construct polling URL (needed for output even if polling fails)
pollURL = [domainRoot '/cgi-bin/webface2.cgi?jobid=' jobId '&wait=' num2str(waitTime)];

% Submit email notification if provided
if ~isempty(email)
    if verbose
        fprintf('  Submitting email notification: %s\n', email);
    end
    
    try
        % URL-encode the email address
        % Use Java URLEncoder for proper encoding (e.g., @ becomes %40)
        if exist('java.net.URLEncoder', 'class') == 8
            emailEncoded = char(java.net.URLEncoder.encode(email, 'UTF-8'));
        else
            % Fallback: manual encoding for common characters
            emailEncoded = email;
            emailEncoded = strrep(emailEncoded, '@', '%40');
            emailEncoded = strrep(emailEncoded, '+', '%2B');
            emailEncoded = strrep(emailEncoded, ' ', '%20');
        end
        
        % Construct email submission URL
        emailURL = [pollURL '&email=' emailEncoded '&submit=Send+email'];
        
        % Submit email via GET request (silent, don't need response)
        webread(emailURL, weboptions('Timeout', timeout));
        
        if verbose
            fprintf('    Email notification submitted successfully\n');
        end
    catch ME
        % Don't fail the whole process if email submission fails
        if verbose
            fprintf('    Warning: Could not submit email notification: %s\n', ME.message);
        end
    end
end

% Poll until CSV is available
maxPollTime = timeout;
pollInterval = waitTime;  % Use the wait time from server
elapsedTime = 0;
csvURL = '';
pollCount = 0;

while elapsedTime < maxPollTime && isempty(csvURL)
    pollCount = pollCount + 1;
    % Print elapsed time every 30 seconds
    if verbose && (mod(elapsedTime, 30) == 0 || pollCount == 1)
        fprintf('    Polling job status (attempt %d, elapsed: %d seconds)...\n', pollCount, elapsedTime);
    end
    
    try
        % Use webread to poll (simpler than matlab.net.http for GET)
        pollHtml = webread(pollURL, weboptions('Timeout', timeout));
        
        % Ensure pollHtml is char (webread may return string)
        if isstring(pollHtml)
            pollHtml = char(pollHtml);
        end
        
        % Check if job is complete and construct CSV URL directly from job ID
        % The CSV URL follows a predictable pattern: /services/DeepLoc-2.1/tmp/{jobId}/results_{jobId}.csv
        % Once the results page is available, the CSV can be downloaded at this URL
        csvURL = '';
        jobComplete = false;
        completionIndicator = '';
        
        % Check for indicators that the job is complete (results page is displayed)
        
        % Check for completion indicators (case-insensitive)
        if contains(pollHtml, 'Download prediction results', 'IgnoreCase', true)
            jobComplete = true;
            completionIndicator = 'Download prediction results';
        elseif contains(pollHtml, 'CSV Summary', 'IgnoreCase', true)
            jobComplete = true;
            completionIndicator = 'CSV Summary';
        elseif contains(pollHtml, 'predicted sequences', 'IgnoreCase', true)
            jobComplete = true;
            completionIndicator = 'predicted sequences';
        elseif contains(pollHtml, 'Finished prediction', 'IgnoreCase', true)
            jobComplete = true;
            completionIndicator = 'Finished prediction';
        end
        
        if jobComplete
            % Job is complete - construct CSV URL directly from job ID
            csvURL = [baseURL '/tmp/' jobId '/results_' jobId '.csv'];
            if verbose
                fprintf('    Detected job completion: Found "%s"\n', completionIndicator);
            end
        end
        
        if ~isempty(csvURL)
            if verbose
                fprintf('    Job completed! CSV URL: %s\n', csvURL);
            end
            break;
        end
        
        % Wait before next poll
        pause(pollInterval);
        elapsedTime = elapsedTime + pollInterval;
        
    catch ME
        EM = sprintf('Error polling job status: %s', ME.message);
        dispEM(EM, false);  % Warning, continue polling
        pause(pollInterval);
        elapsedTime = elapsedTime + pollInterval;
    end
end

if isempty(csvURL)
    EM = sprintf('Job did not complete within timeout (%d seconds). Job ID: %s', timeout, jobId);
    if verbose
        fprintf('\n');
        fprintf('  Job ID: %s\n', jobId);
        fprintf('  Number of polling attempts: %d\n', pollCount);
        fprintf('  You can manually check the job status by visiting:\n');
        fprintf('    %s\n', pollURL);
    end
    dispEM(EM, true);
end

% ============================================================================
% STEP 4: Download CSV file
% ============================================================================
if verbose
    fprintf('  Step 4: Downloading CSV file: %s\n', csvURL);
end

try
    % Use webread to download CSV
    % Note: webread may automatically parse CSV into a table, so we need to handle that
    csvContent = webread(csvURL, weboptions('Timeout', timeout));
    
    % webread may automatically parse CSV into a table, so handle that case
    if istable(csvContent)
        % Convert table back to CSV string format
        tempFile = [tempname '.csv'];
        writetable(csvContent, tempFile);
        csvContent = fileread(tempFile);
        delete(tempFile);
    elseif isstring(csvContent)
        csvContent = char(csvContent);
    end
    
    if verbose
        fprintf('    CSV file downloaded (%d characters)\n', length(csvContent));
    end
catch ME
    EM = sprintf('Error downloading CSV file: %s', ME.message);
    dispEM(EM, true);
end

% Generate output filename
[~, fastaName, ~] = fileparts(fastaFilePath);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csvFilename = sprintf('deeploc_%s_%s.csv', fastaName, timestamp);
csvFilePath = fullfile(outputDir, csvFilename);

% Save CSV file
try
    fid = fopen(csvFilePath, 'w');
    if fid == -1
        EM = sprintf('Cannot create CSV file: %s', csvFilePath);
        dispEM(EM, true);
    end
    fprintf(fid, '%s', csvContent);
    fclose(fid);
catch ME
    EM = sprintf('Error saving CSV file: %s', ME.message);
    dispEM(EM, true);
end

if verbose
    fprintf('Success! CSV file saved to: %s\n', csvFilePath);
end

% Return pollURL as second output
% pollURL is already defined from the polling loop

end