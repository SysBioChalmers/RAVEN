function [gffFile, faaFile] = downloadGenomeData(accession, outputDir, verbose)
% downloadGenomeData  Download genome annotation files from NCBI FTP.
%
% Retrieves the GFF3 annotation and protein FASTA (.faa) files for a given
% NCBI genome assembly accession using the NCBI E-utilities API to resolve
% the FTP path automatically. The files are saved locally and can be used
% directly by getGeneData to build a gene-mapping table.
%
% Parameters
% ----------
% accession : char
%     NCBI genome assembly accession, e.g. 'GCF_000002595.2'.
%     Both RefSeq (GCF_) and GenBank (GCA_) prefixes are accepted.
% outputDir : char, optional
%     Directory where downloaded files are saved (default: current working
%     directory).
% verbose : logical, optional
%     Print download progress to the command window (default: true).
%
% Returns
% -------
% gffFile : char
%     Full path to the downloaded GFF3 annotation file.
% faaFile : char
%     Full path to the downloaded protein FASTA file.
%
% Examples
% --------
%     % Download to current directory
%     [gff, faa] = downloadGenomeData('GCF_000002595.2');
%
%     % Save to a specific directory
%     [gff, faa] = downloadGenomeData('GCF_000002595.2', 'data/');
%
%     % Suppress progress messages
%     [gff, faa] = downloadGenomeData('GCF_000002595.2', 'data/', false);
%
% Notes
% -----
% Requires an active internet connection. Files already present in
% outputDir are not re-fetched; delete them manually to force a refresh.

    if nargin < 1 || isempty(accession)
    error('downloadGenomeData:missingInput', ...
        ['accession is required. ', ...
         'Provide an NCBI assembly accession, e.g. ''GCF_000002595.2''.']);
    end
    accession = char(accession);
    
    if ~any(strncmpi(accession, {'GCF_', 'GCA_'}, 4))
        error('downloadGenomeData:invalidAccession', ...
        ['Accession must start with ''GCF_'' (RefSeq) or ''GCA_'' (GenBank). ', ...
         'Received: ''%s''.'], accession);
    end

    if nargin < 2 || isempty(outputDir)
        downloadDir = pwd;
    else
        downloadDir = char(outputDir);
    end

    if nargin < 3
        verbose = true;
    end
    
    % Use calls with conditional logging instead of fprintf
	logMsg = @(varargin) verbose && fprintf(varargin{:});

	logMsg('Looking up NCBI FTP path for accession: %s\n', accession);
    if ~exist(downloadDir, 'dir')
        mkdir(downloadDir);
    end
    
    % --- E-utilities: search for the assembly ---
    baseAcc = regexprep(accession, '\.\d+$', '');  % strip version for search
    searchUrl = sprintf( ...
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=%s&retmode=json', ...
        baseAcc);
    try
        searchData = webread(searchUrl);
    catch ME
        error('getGeneData:download', 'NCBI E-search failed: %s', ME.message);
    end
    
    if isempty(searchData.esearchresult.idlist)
        error('getGeneData:download', 'No assembly found for accession: %s', accession);
    end
    uid = searchData.esearchresult.idlist{1};
    
    % --- E-utilities: fetch assembly summary to get FTP path ---
    summaryUrl = sprintf( ...
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=%s&retmode=json', uid);
    try
        summaryData = webread(summaryUrl);
    catch ME
        error('getGeneData:download', 'NCBI E-summary failed: %s', ME.message);
    end
    
    docsum  = summaryData.result.(sprintf('x%s', uid));
    
    % GCF_ accessions are RefSeq assemblies; GCA_ are GenBank assemblies.
    % Use the matching FTP path first, fall back to the other if empty.
    if strncmpi(accession, 'GCF', 3)
        ftpBase    = docsum.ftppath_refseq;
        ftpFallback = docsum.ftppath_genbank;
    else  % GCA_
        ftpBase    = docsum.ftppath_genbank;
        ftpFallback = docsum.ftppath_refseq;
    end
    
    if isempty(ftpBase)
        warning('getGeneData:ftpFallback', ...
            'Primary FTP path not found for %s, trying fallback.', accession);
        ftpBase = ftpFallback;
    end
    if isempty(ftpBase)
        error('getGeneData:download', 'No FTP path found for accession: %s', accession);
    end
    
    % Build base names and URLs for GFF3 and protein FASTA (NCBI naming:
    % <folder_name>_genomic.gff.gz and <folder_name>_protein.faa.gz)
    pathParts    = strsplit(ftpBase, '/');
    folderName   = pathParts{end};
    % Construct HTTPS URLs (use https mirror of NCBI FTP)
    gffGzUrl     = strrep(sprintf('%s/%s_genomic.gff.gz',  ftpBase, folderName), 'ftp://', 'https://');
    proteinGzUrl = strrep(sprintf('%s/%s_protein.faa.gz',  ftpBase, folderName), 'ftp://', 'https://');
    % Local paths
    genomePath  = downloadAndExtract(gffGzUrl,     downloadDir, [folderName '_genomic.gff.gz'], verbose);
    proteinPath = downloadAndExtract(proteinGzUrl, downloadDir, [folderName '_protein.faa.gz'], verbose);

    gffFile  = genomePath;
	faaFile  = proteinPath;
end

%--------------------------------------------------------------------------
% Helper functions

function outFile = downloadAndExtract(url, destDir, gzName, verbose)
% Downloads a .gz file from url into destDir, extracts it, and returns
% the path to the extracted file. Skips download if already extracted.
    gzFile  = fullfile(destDir, gzName);
    outFile = fullfile(destDir, regexprep(gzName, '\.gz$', ''));

    if isfile(outFile)
        if verbose; fprintf('%s already present, skipping download.\n', gzName); end
        return
    end

    if verbose; fprintf('Downloading: %s\n', url); end
    try
        websave(gzFile, url);
    catch ME
        error('downloadGenomeData:download', 'Download failed for %s: %s', url, ME.message);
    end
    try
        gunzip(gzFile, destDir);
        delete(gzFile);
    catch ME
        error('downloadGenomeData:extract', 'Extraction failed for %s: %s', gzName, ME.message);
    end
end