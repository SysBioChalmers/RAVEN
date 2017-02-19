function exitFlag=downloadKEGG(keggPath)
% downloadKEGG
%   Retrieves files used by getKEGGModelForOrganism from the KEGG FTP
%
%   keggPath    the directory to store the files in
%
%   exitFlag    1 if everything worked, -1 if something went wrong
%
%   Usage: exitFlag=downloadKEGG(keggPath)
%
%   Rasmus Agren, 2014-01-09
%

exitFlag=-1;
%Check if all necessary files exist in the KEGG path. Otherwise download
%them from KEGG
necessaryFiles={'compound';'compound.inchi';'genes.pep';'ko';'reaction';'reaction.lst';'taxonomy';'reaction_mapformula.lst'};
filePathways={'/pub/kegg/ligand/compound/';'/pub/kegg/ligand/compound/';'/pub/kegg/genes/fasta/';'/pub/kegg/genes/';'/pub/kegg/ligand/reaction/';'/pub/kegg/ligand/reaction/';'/pub/kegg/genes/';'pub/kegg/ligand/reaction/'};
fileExists=false(numel(necessaryFiles),1);
for i=1:numel(necessaryFiles)
    if exist(fullfile(keggPath,necessaryFiles{i}),'file')
        fileExists(i)=true;
    end
end

%Download missing files
if any(~fileExists)
    retrieveFiles=input('One or more files are missing from the local database. Do you want to connect to KEGG and retrieve them? (Y/N)','s');
    if strcmpi(retrieveFiles,'y')
        fprintf('\nWARNING: This might take a very long time. The database is about 2 GB.\n');
        getFiles=find(~fileExists);
        fprintf('Attempting to connect to the KEGG database and retrieve the files. The FTP client in Matlab has some known issues. If you do not see that it is sucessfully downloading the files you could get them manually from KEGG. Files to get and place in keggPath:\n');
        for i=1:numel(getFiles)
            fprintf(['ftp.genome.jp' filePathways{getFiles(i)} necessaryFiles{getFiles(i)} '\n']);
        end

        f = ftp('ftp.genome.jp');
        pasv(f);
        for i=numel(getFiles)
            cd(f,filePathways{getFiles(i)});
            fprintf(['Retrieving ' filePathways{getFiles(i)} necessaryFiles{getFiles(i)} '\n']);
            mget(f,necessaryFiles{getFiles(i)},keggPath);
        end
        close(ftp);
        exitFlag=1;
    end
end
end
