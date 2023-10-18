function GSS = getWoLFScores(inputFile, kingdom)
% getWoLFScores
%   Call WoLF PSort to predict the sub-cellular localization of proteins.
%   The output can be used as input to predictLocalization. This function
%   is currently only available for Linux and requires Perl to be
%   installed. If one wants to use another predictor, see parseScores. The
%   function normalizes the scores so that the best score for each gene is
%   1.0.
%
%   Input:
%   inputFile	a FASTA file with protein sequences
%   kingdom     the kingdom of the organism, 'animal', 'fungi' or 'plant'
%
%   Output:
%   GSS         a gene scoring structure to be used in predictLocalization
%
%   Usage: GSS = getWoLFScores(inputFile, kingdom)

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

%Do the prediction
[~, output]=unix(['perl "' ravenPath '/software/WoLFPSORT/bin/runWolfPsortSummary" ' kingdom ' < ' inputFile]);

%Save output and call the general parser
fprintf(fid,output);
fclose(fid);
GSS=parseScores(outFile,'wolf');

%Clean up
delete(outFile);
end
