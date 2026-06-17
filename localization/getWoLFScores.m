function GSS = getWoLFScores(inputFile, kingdom)
% getWoLFScores  Predict protein sub-cellular localization with WoLF PSORT.
%
% The output can be used as input to predictLocalization. This function is
% currently only available for Linux and requires Perl to be installed. If
% one wants to use another predictor, see parseScores. The function
% normalizes the scores so that the best score for each gene is 1.0.
%
% Parameters
% ----------
% inputFile : char
%     a FASTA file with protein sequences.
% kingdom : char
%     the kingdom of the organism, 'animal', 'fungi' or 'plant'.
%
% Returns
% -------
% GSS : struct
%     a gene scoring structure to be used in predictLocalization.
%
% Examples
% --------
%     GSS = getWoLFScores(inputFile, kingdom);
%
% See also
% --------
% parseScores, predictLocalization

if ~isfile(inputFile)
    error('FASTA file %s cannot be found',string(inputFile));
end

kingdom=char(kingdom);
if ~any(strcmp(kingdom,{'animal','fungi','plant'}))
    EM='Allowed kingdoms are "animal", "fungi", and "plant"';
    dispEM(EM);
end

if ispc || ismac
    EM='This function currently runs only on Linux. Use parseScores if you want to use another predictor';
    dispEM(EM);
end

%Get the directory for RAVEN Toolbox
ravenPath=findRAVENroot();

%Temporary output name
outFile=tempname;
fid=fopen(outFile,'w');

%Fetch the WoLFPSORT bundle on demand if it is not already present
if ~exist(fullfile(ravenPath,'software','WoLFPSORT','bin','runWolfPsortSummary'),'file')
    downloadRavenBinaries({'WoLFPSORT'});
end

%Do the prediction
[~, output]=unix(['perl "' ravenPath '/software/WoLFPSORT/bin/runWolfPsortSummary" ' kingdom ' < ' inputFile]);

%Save output and call the general parser
fprintf(fid,output);
fclose(fid);
GSS=parseScores(outFile,'wolf');

%Clean up
delete(outFile);
end
