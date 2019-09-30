function geneScoreStruct=getWoLFScores(inputFile,kingdom)
% getWoLFScores
%   Calls WoLF PSort to predict the sub-cellular localization of proteins.
%   The output can be used as input to predictLocalization. This function
%   is currently only available for Linux and requires PERL to be
%   installed. If you want to use another predictor, see parseScores.
%   The function normalizes the scores so that the best score for each gene
%   is 1.0.
%
%   Input:
%   inputFile               a FASTA file with protein sequences
%   kingdom                 the kingdom of the organism, 'animal', 'fungi'
%                           or 'plant'
%
%   Output:
%   geneScoreStructure      a structure to be used in predictLocalization
%
%   Usage: geneScoreStruct=getWoLFScores(inputFile,kingdom)
%
%   Simonas Marcisauskas, 2019-08-17
%

if ~(exist(inputFile,'file')==2)
    error('FASTA file %s cannot be found',string(inputFile));
end

if ~strcmp(kingdom,'animal') && ~strcmp(kingdom,'fungi') && ~strcmp(kingdom,'plant')
    EM='Allowed kingdoms are "animal", "fungi", and "plant"';
    dispEM(EM);
end

if ispc==true
    EM='This function currently runs only on Linux. Use parseScores if you want to use another predictor';
    dispEM(EM);
end

%Get the directory for RAVEN Toolbox. This may not be the easiest or best
%way to do this
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));

%Temporary output name
outFile=tempname;
fid=fopen(outFile,'w');

%Do the prediction
[~, output]=unix(['perl "' ravenPath '/software/WoLFPSORT/bin/runWolfPsortSummary" ' kingdom ' < ' inputFile]);

%Save output and call the general parser
fprintf(fid,output);
fclose(fid);
geneScoreStruct=parseScores(outFile,'wolf');

%Clean up
delete(outFile);
end
