# DeepLoc Integration for RAVEN Toolbox

## Overview

This module provides integration with DeepLoc and DeepLocPro webservers for protein subcellular localization prediction. The results are automatically formatted for use with RAVEN's `predictLocalization()` function.

## Workflow

```
FASTA file → deeplocRun.m / deeplocproRun.m → CSV file → parseScores('deeploc') → GSS → predictLocalization.m
```

## Quick Start

### DeepLoc (Eukaryotic proteins)

```matlab
% Step 1: Run DeepLoc and get CSV results
[csvFile, pollURL] = deeplocRun('proteins.fasta');

% Step 2: Parse CSV to get GSS structure
GSS = parseScores(csvFile, 'deeploc');

% Step 3: Use GSS with predictLocalization
[model, geneLocalization, transportStruct, scores, removedRxns] = ...
    predictLocalization(model, GSS, 'Cytoplasm', 0.5, 15, false);
```

### DeepLocPro (Prokaryotic proteins)

```matlab
% Step 1: Run DeepLocPro and get CSV results
[csvFile, pollURL] = deeplocproRun('proteins.fasta', 'OrganismGroup', 'Gram negative');

% Step 2: Parse CSV to get GSS structure
GSS = parseScores(csvFile, 'deeploc');

% Step 3: Use GSS with predictLocalization
[model, geneLocalization, transportStruct, scores, removedRxns] = ...
    predictLocalization(model, GSS, 'Cytoplasm', 0.5, 15, false);
```

## Functions

### `deeplocRun.m`

Main function that submits FASTA sequences to DeepLoc webserver (versions 1.0, 2.0, or 2.1) and downloads CSV results.

**Usage:**
```matlab
[csvFile, pollURL] = deeplocRun('proteins.fasta');
[csvFile, pollURL] = deeplocRun('proteins.fasta', 'Version', '2.0');
[csvFile, pollURL] = deeplocRun('proteins.fasta', 'Version', '1.0', 'Encoding', 'blosum62');
[csvFile, pollURL] = deeplocRun('proteins.fasta', 'Mode', 'fast', 'Figures', true, 'Email', 'user@example.com');
```

**Options:**
- `'Version'` - DeepLoc version: '2.1' (default), '2.0', or '1.0'
- `'Mode'` - prediction mode: 'high-quality' (default) or 'fast' (for DeepLoc 2.0/2.1 only)
- `'Figures'` - true/false (default: false) - request figures from DeepLoc (for DeepLoc 2.0/2.1 only)
- `'Encoding'` - protein encoding: 'profiles' (default) or 'blosum62' (for DeepLoc 1.0 only)
- `'Email'` - email address for notification when job completes (optional)
- `'OutputDir'` - directory to save CSV file (default: same directory as FASTA file)
- `'Timeout'` - HTTP request timeout in seconds (default: 900)
- `'Verbose'` - print progress messages (default: true)

**Output:**
- `csvFilePath` - full path to downloaded CSV file with DeepLoc results
- `pollURL` - URL to check job status on DeepLoc webserver

### `deeplocproRun.m`

Main function that submits prokaryotic protein sequences to DeepLocPro-1.0 webserver and downloads CSV results.

**Usage:**
```matlab
[csvFile, pollURL] = deeplocproRun('proteins.fasta');
[csvFile, pollURL] = deeplocproRun('proteins.fasta', 'OrganismGroup', 'Gram negative');
[csvFile, pollURL] = deeplocproRun('proteins.fasta', 'OrganismGroup', 'Archaea', 'Figures', true, 'Email', 'user@example.com');
```

**Options:**
- `'OrganismGroup'` - organism group: 'Any' (default), 'Archaea', 'Gram negative', or 'Gram positive'
- `'Figures'` - true/false (default: false) - request figures from DeepLocPro
- `'Email'` - email address for notification when job completes (optional)
- `'OutputDir'` - directory to save CSV file (default: same directory as FASTA file)
- `'Timeout'` - HTTP request timeout in seconds (default: 900)
- `'Verbose'` - print progress messages (default: true)

**Output:**
- `csvFilePath` - full path to downloaded CSV file with DeepLocPro results
- `pollURL` - URL to check job status on DeepLocPro webserver

## Requirements

- Internet connection (for DeepLoc/DeepLocPro webservers)
- MATLAB R2023b or later (uses `matlab.net.http` package)
- RAVEN Toolbox installed

## How It Works

1. **Submission**: Reads FASTA file and submits it to DeepLoc/DeepLocPro webserver via HTTP POST
2. **Job Tracking**: Captures job ID from server redirect response (302 Found)
3. **Polling**: Polls job status until completion is detected
4. **Download**: Extracts CSV URL from JSON endpoint (DeepLocPro) or constructs from job ID (DeepLoc) and downloads results
5. **Output**: Saves CSV file to specified output directory

## Debugging

For debugging purposes, you can use `deeplocRun_Debug.m` or `deeplocproRun_Debug.m` which break down the functions into step-by-step sections that can be run interactively.

## Notes

- The CSV format is compatible with `parseScores.m` using the 'deeploc' option
- This module does not modify any core RAVEN functions
- DeepLoc results page uses AngularJS templates; CSV URL is constructed from job ID
- DeepLocPro CSV URL includes a date-time suffix and is extracted from JSON endpoint
- Maximum sequence limits: 500 sequences (DeepLoc 2.0/2.1, DeepLocPro), 50 sequences (DeepLoc 1.0 with Profiles encoding), 500 sequences (DeepLoc 1.0 with BLOSUM62 encoding)
