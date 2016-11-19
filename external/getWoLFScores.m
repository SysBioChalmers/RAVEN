function geneScoreStruct=getWoLFScores(inputFile,kingdom)
% getWoLFScores
%   Calls WoLF PSort to predict the sub-cellular localization of proteins.
%   The output can be used as input to predictLocalization. This function
%   is currently only available for Linux and requires PERL to be
%   installed. If you want to use another predictor, see parseScores.
%
%   inputFile	a FASTA file with protein sequences
%   kingdom     the kingdom of the organism, 'animal', 'fungi' or 'plant'.
%
%   The function normalizes the scores so that the best score for each gene
%   is 1.0.
%
%   geneScoreStructure  a structure to be used in predictLocalization
%
%   Usage: geneScoreStruct=getWoLFScores(inputFile,kingdom)
%
%   Rasmus Agren, 2012-03-27

if ~strcmp(kingdom,'animal') && ~strcmp(kingdom,'fungi') && ~strcmp(kingdom,'plant')
    dispEM('Allowed kingdoms are "animal", "fungi", and "plant"'); 
end

if ispc==true
    dispEM('This function currently runs only on Linux. Use parseScores if you want to use another predictor'); 
end

%Get the directory for RAVEN Toolbox. This may not be the easiest or best
%way to do this
[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(ST(I).file));
% Adding escape characters, if some parent folders contain spaces or
% exclamation marks (for Unix systems). For Windows, all the parent folders
% are just put between the double quotation brackets
if isunix
    ravenPath = regexprep(ravenPath,'\ ','\\ ');
    ravenPath = regexprep(ravenPath,'\!','\\!');
elseif ispc
    for i=1:(length(strfind(ravenPath,'\')))
        if i==1
            ravenPath = regexprep(ravenPath,'\\','\\"',i);
        elseif i==length(strfind(ravenPath,'\'))
            ravenPath = regexprep(ravenPath,'\\','"\\',i);    
        else
            ravenPath = regexprep(ravenPath,'\\','"\\"',i);
        end
    end
end

%Temporary output name
outFile=tempname;
fid=fopen(outFile,'w');

%Do the prediction
[crap output]=unix(['perl ' ravenPath '/software/WoLFPSORT_package_v0.2/bin/runWolfPsortSummary ' kingdom ' < ' inputFile]);

%Save output and call the general parser
fprintf(fid,output);
fclose(fid);
geneScoreStruct=parseScores(outFile,'wolf');

%Clean up
delete(outFile);
