function [csvFilePath, pollURL] = deeplocproRun(fastaFilePath, varargin)
% deeplocproRun
%   Submits prokaryotic protein sequences from a FASTA file to DeepLocPro
%   webserver and downloads the resulting CSV file with localization predictions.
%
% Input:
%   fastaFilePath     path to FASTA file containing protein sequences
%   'OrganismGroup'   organism group: 'Any' (default), 'Archaea', 'Gram negative',
%                     or 'Gram positive' (optional)
%   'Figures'         true if figures should be requested from DeepLocPro
%                     (optional, default false)
%   'Email'           email address for notification when job completes
%                     (optional, default: no email)
%   'OutputDir'       directory to save CSV file (optional, default: same
%                     directory as FASTA file)
%   'Timeout'         HTTP request timeout in seconds (optional, default 900)
%   'Verbose'         true if progress messages should be printed (optional,
%                     default true)
%
% Output:
%   csvFilePath       full path to downloaded CSV file with DeepLocPro results
%   pollURL           URL to check job status on DeepLocPro webserver
%
% Usage:
%   [csvFile, pollURL] = deeplocproRun('proteins.fasta');
%   [csvFile, pollURL] = deeplocproRun('proteins.fasta', 'OrganismGroup', 'Gram negative');
%   [csvFile, pollURL] = deeplocproRun('proteins.fasta', 'OrganismGroup', 'Archaea', ...
%       'Figures', true, 'Email', 'user@example.com');
%
% NOTE: This function requires internet connection to access DeepLocPro-1.0
% webserver. The CSV file format must be compatible with parseScores.m
% using the 'deeploc' option (or 'deeplocpro' if implemented).
%
% The workflow is:
%   1. deeplocproRun() → CSV file
%   2. parseScores(CSV, 'deeploc') → GSS structure
%   3. predictLocalization(model, GSS, ...) → compartmentalized model

% Valid organism groups
validOrganismGroups = {'Any', 'Archaea', 'Gram negative', 'Gram positive'};

% Parse optional parameters
p = inputParser;
addParameter(p, 'OrganismGroup', 'Any', @(x) ischar(x) || isstring(x));
addParameter(p, 'Figures', false, @islogical);
addParameter(p, 'Email', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x) || isempty(x));  % Will be set to FASTA directory if empty
addParameter(p, 'Timeout', 900, @isnumeric);
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});

organismGroup = char(p.Results.OrganismGroup);
figures = p.Results.Figures;
email = char(p.Results.Email);
outputDir = char(p.Results.OutputDir);
timeout = p.Results.Timeout;
verbose = p.Results.Verbose;

% Validate organism group
if ~ismember(organismGroup, validOrganismGroups)
    EM = sprintf('OrganismGroup ''%s'' is not supported. Valid groups are: %s', ...
        organismGroup, strjoin(validOrganismGroups, ', '));
    dispEM(EM, true);
end

% Map organism group to DeepLocPro values (lowercase with underscores)
if strcmpi(organismGroup, 'Any')
    groupValue = 'any';
elseif strcmpi(organismGroup, 'Archaea')
    groupValue = 'archaea';
elseif strcmpi(organismGroup, 'Gram negative')
    groupValue = 'gram_negative';
elseif strcmpi(organismGroup, 'Gram positive')
    groupValue = 'gram_positive';
else
    EM = sprintf('Invalid organism group: %s', organismGroup);
    dispEM(EM, true);
end

% Map figures parameter to format
if figures
    formatValue = 'long';
else
    formatValue = 'short';
end

% Construct baseURL for DeepLocPro
baseURL = 'https://services.healthtech.dtu.dk/services/DeepLocPro-1.0';

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

% Try to read FASTA file to validate format
try
    % Try using fastaread if Bioinformatics Toolbox is available
    if exist('fastaread', 'file') == 2
        sequences = fastaread(fastaFilePath);
        if isempty(sequences)
            EM = 'FASTA file appears to be empty';
            dispEM(EM, true);
        end
        numSequences = numel(sequences);
        
        % Validate sequence count (max 500 for DeepLocPro)
        if numSequences > 500
            EM = sprintf('DeepLocPro accepts a maximum of 500 sequences. Your file contains %d sequences.', numSequences);
            dispEM(EM, true);
        end
        
        % Validate sequence lengths (10-6000 amino acids)
        shortSequences = [];
        longSequences = [];
        for i = 1:numSequences
            seqLen = length(sequences(i).Sequence);
            if seqLen < 10
                shortSequences = [shortSequences, i];
            elseif seqLen > 6000
                longSequences = [longSequences, i];
            end
        end
        
        if ~isempty(shortSequences)
            EM = sprintf('Warning: %d sequence(s) shorter than 10 amino acids detected. DeepLocPro predictions may be inaccurate for these sequences.', length(shortSequences));
            dispEM(EM, false);
        end
        if ~isempty(longSequences)
            EM = sprintf('Warning: %d sequence(s) longer than 6000 amino acids detected. DeepLocPro predictions may be inaccurate for these sequences.', length(longSequences));
            dispEM(EM, false);
        end
        
        if verbose
            fprintf('  Found %d sequence(s)\n', numSequences);
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
        if numSequences > 500
            EM = sprintf('DeepLocPro accepts a maximum of 500 sequences. Your file appears to contain %d sequences.', numSequences);
            dispEM(EM, true);
        end
        if verbose
            fprintf('  Found %d sequence(s)\n', numSequences);
        end
    end
catch ME
    EM = sprintf('Error reading FASTA file: %s', ME.message);
    dispEM(EM, true);
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
    fprintf('Submitting to DeepLocPro-1.0...\n');
    fprintf('  Organism group: %s (group=%s)\n', organismGroup, groupValue);
    fprintf('  Format: %s\n', formatValue);
end

% ============================================================================
% STEP 1: Manually construct multipart/form-data body
% ============================================================================
if verbose
    fprintf('  Step 1: Preparing multipart form data...\n');
end

% Construct configfile path
configfile = '/var/www/services/services/DeepLocPro-1.0/webface.cf';

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

% Part 2: fasta (pasted sequences - read from file)
% Read FASTA file content as text for pasting
fid = fopen(fastaFilePath, 'r');
if fid == -1
    EM = sprintf('Cannot open FASTA file for reading: %s', fastaFilePath);
    dispEM(EM, true);
end
fastaTextContent = fread(fid, '*char')';
fclose(fid);
% Ensure fastaTextContent is a row vector
if size(fastaTextContent, 1) > 1
    fastaTextContent = fastaTextContent';
end

% Normalize line endings in FASTA content to \r\n (if needed)
% The server expects \r\n line endings
fastaTextContent = regexprep(fastaTextContent, '\r?\n', '\r\n');

% Remove ALL trailing whitespace/newlines from FASTA content
fastaTextContent = regexprep(fastaTextContent, '[\r\n\s]+$', '');

bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="fasta"\r\n');
bodyParts{end+1} = sprintf('\r\n');
% Insert FASTA text content
bodyParts{end+1} = fastaTextContent;
% Add \r\n after FASTA content (before next boundary) - this must be a separate bodyPart
bodyParts{end+1} = sprintf('\r\n');

% Part 3: uploadfile (empty file upload)
bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="uploadfile"; filename=""\r\n');
bodyParts{end+1} = sprintf('Content-Type: application/octet-stream\r\n');
bodyParts{end+1} = sprintf('\r\n');  % Empty field (just one blank line)
bodyParts{end+1} = sprintf('\r\n');  % Second blank line before next boundary

% Part 4: group (organism group)
bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="group"\r\n');
bodyParts{end+1} = sprintf('\r\n');
bodyParts{end+1} = sprintf('%s\r\n', groupValue);

% Part 5: format
bodyParts{end+1} = sprintf('--%s\r\n', boundary);
bodyParts{end+1} = sprintf('Content-Disposition: form-data; name="format"\r\n');
bodyParts{end+1} = sprintf('\r\n');
bodyParts{end+1} = sprintf('%s\r\n', formatValue);

% Closing boundary
bodyParts{end+1} = sprintf('--%s--\r\n', boundary);

% Convert all parts to uint8 and concatenate
bodyChar = [bodyParts{:}];
bodyUint8 = uint8(bodyChar);
% Ensure bodyUint8 is a row vector
if size(bodyUint8, 1) > 1
    bodyUint8 = bodyUint8';
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
    EM = sprintf('Error submitting POST request to DeepLocPro: %s', ME.message);
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
        
        % Check if job is complete and extract CSV URL from JSON endpoint
        csvURL = '';
        jobComplete = false;
        completionIndicator = '';
        
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
            % Job is complete - extract CSV URL from JSON endpoint
            % Fetch JSON data from known endpoints to get the CSV file URL
            jsonEndpoints = {
                [baseURL '/tmp/' jobId '/results.json'],
                [baseURL '/tmp/' jobId '/data.json'],
                [domainRoot '/services/DeepLocPro-1.0/tmp/' jobId '/results.json'],
                [domainRoot '/services/DeepLocPro-1.0/tmp/' jobId '/data.json']
            };
            
            for e = 1:length(jsonEndpoints)
                try
                    jsonData = webread(jsonEndpoints{e}, weboptions('Timeout', 10));
                    if isstruct(jsonData) && isfield(jsonData, 'csv_file')
                        csvURL = jsonData.csv_file;
                        if ~startsWith(csvURL, 'http')
                            if startsWith(csvURL, '/')
                                csvURL = ['https://services.healthtech.dtu.dk' csvURL];
                            else
                                csvURL = [baseURL '/' csvURL];
                            end
                        end
                        if verbose
                            fprintf('    Found CSV URL from JSON endpoint %d: %s\n', e, csvURL);
                        end
                        break;
                    end
                catch
                    % Endpoint doesn't exist or failed, continue to next
                end
            end
            
            if ~isempty(csvURL)
                if verbose
                    fprintf('    Detected job completion: Found "%s"\n', completionIndicator);
                end
            else
                % If we can't find the CSV URL
                if verbose
                    fprintf('    Detected job completion: Found "%s" but could not extract CSV URL from JSON endpoints\n', completionIndicator);
                end
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
    EM = sprintf('DeepLocPro job did not complete within timeout (%d seconds). Job ID: %s', timeout, jobId);
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
    csvContent = webread(csvURL, weboptions('Timeout', timeout));
    
    % webread may automatically parse CSV files into tables
    % We need to handle this and convert back to text if needed
    if istable(csvContent)
        % If webread parsed it as a table, write to temp file and read back as text
        tempFile = [tempname '.csv'];
        writetable(csvContent, tempFile);
        csvContent = fileread(tempFile);
        delete(tempFile);
    end
    
    % Ensure csvContent is char (webread may return string)
    if isstring(csvContent)
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
csvFilename = sprintf('deeplocpro_%s_%s.csv', fastaName, timestamp);
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

end

