function constructMultiFasta(model,sourceFile,outputDir)
% constructMultiFasta
%   Saves one file in FASTA format for each reaction in the model that has genes
%
%   Input:
%   model         a model structure
%   sourceFile    a file with sequences in FASTA format
%   outputDir     the directory to save the resulting FASTA files in
%
%   The source file is assumed to have the format '>gene identifier
%   additional info'. Only the gene identifier is used for matching. This is
%   to be compatible with the rest of the code that retrieves information
%   from KEGG.
%
%   Usage: constructMultiFasta(model,sourceFile,outputDir)

sourceFile=char(sourceFile);
outputDir=char(outputDir);
if ~(exist(sourceFile,'file')==2)
    error('FASTA file %s cannot be found',string(sourceFile));
end

fprintf('Scanning the source multi-FASTA file... ');
%Open the source file
fid = fopen(sourceFile, 'r');

%Since the list of genes will be accessed many times I use a Java hashtable
hTable = java.util.Hashtable;

for i=1:numel(model.genes)
    hTable.put(model.genes{i}, i);
end

%First go through the source file and save the position (in bytes) of the
%start of each gene
elementPositions=zeros(5000000,1,'int64'); %Make room for 5000000 elements
totalElements=0;
whereAmI=0; %Keeps track of where in the file we are
while 1
    %Read 10 mb at a time
    str=fread(fid,10000000,'int8');
    
    %Find any '>' which indicates a new label in FASTA format
    newPosts=find(str==62); %62 is '>'
    
    elementPositions(totalElements+1:totalElements+numel(newPosts))=whereAmI+newPosts;
    
    totalElements=totalElements+numel(newPosts);
    
    whereAmI=whereAmI+10000000;
    
    if feof(fid)
        break;
    end
end
elementPositions=elementPositions(1:totalElements);
fprintf('COMPLETE\n');

fprintf(['NOTICE: If Matlab is freezing and does not provide any output in 30 minutes, consider increasing Java Heap Memory\n', ...
    'in MATLAB settings and start over with the new session\n']);
fprintf('Mapping genes to the multi-FASTA source file...   0%% complete');
%Now loop through the file to see which genes are present in the gene list
%and save their position IN elementPositions! This is to enable a easy way
%to get the distance to the following element
genePositions=zeros(numel(model.genes),1);
for i=1:numel(elementPositions)
    fseek(fid,elementPositions(i),-1);
    str=fread(fid,[1 30],'*char'); %Assumes that no ID is longer than 20 characters
    delim=find(str==32 | str==10,1,'first'); %Space or line feed
    
    geneIdentifier=str(1:delim-1);
    
    %This should never happen, but just to prevent errors. Could be that
    %'>' is a part of some gene information. An alternative would be to
    %check that the indexes follows a line feed
    if isempty(geneIdentifier)
        continue;
    end
    
    %If not found it means that the id was too long
    if isempty(delim)
        EM='Too long gene identifier, increase read length';
        dispEM(EM);
    end
    
    %See if the gene was found
    id=hTable.get(geneIdentifier);
    
    if any(id)
        if genePositions(id)==0
            genePositions(id)=i;
        end
    end
    %Print the progress
    if rem(i-1,350000) == 0
        progress=num2str(floor(100*i/numel(elementPositions)));
        progress=pad(progress,3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');   

fprintf('Generating the KEGG Orthology specific multi-FASTA files...   0%% complete');
%Loop through the reactions and print the corresponding sequences
for i=1:numel(model.rxns)

    %Do not overwrite existing files
    if ~exist(fullfile(outputDir,[model.rxns{i} '.fa']), 'file')
        
        %Get the positions in elementPositions for the involved genes
        genesUsed=model.rxnGeneMat(i,:);
        
        %Open a file for this reaction. This saves empty files for KOs
        %without genes
        rxnfid=fopen(fullfile(outputDir,[model.rxns{i} '.fa']),'w');
        
        if any(genesUsed)
            positions=genePositions(genesUsed~=0);
            
            %It could be that some genes were not found. Delete those
            %elements
            positions(positions==0)=[];
            
            %Print each sequence
            for j=1:numel(positions)
                fseek(fid,elementPositions(positions(j)),-1);
                %Should check that it ends with a gene!!!! Check for eof
                if positions(j)<numel(elementPositions)
                    str=fread(fid,[1 double(elementPositions(positions(j)+1))-double(elementPositions(positions(j)))-1],'*char');
                    
                    %If the string does not end with a line feed character
                    if str(end)~=10
                        str=[str fread(fid,[1 double(elementPositions(positions(j)+2))-double(elementPositions(positions(j)+1))],'*char')];
                        
                        %This is if we still have not found a new gene.
                        %Maximal unluck. This whole check should be done
                        %when elementPositions are calculated!
                        if str(end)~=10
                            %Skip this gene
                            continue;
                        end
                    end
                else
                    str=fread(fid,[1 inf],'*char');
                end
                fwrite(rxnfid,['>' str]);
            end
        end
        fclose(rxnfid);
    end
    %Print the progress
    if rem(i-1,50) == 0
        progress=num2str(floor(100*i/numel(model.rxns)));
        progress=pad(progress,3,'left');
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%s%% complete',progress);
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bCOMPLETE\n');

%Close the source file
fclose(fid);
end

function files=listFiles(directory)
%Supporter function to list the files in a directory and return them as a
%cell array
temp=dir(directory);
files=cell(numel(temp),1);
for i=1:numel(temp)
    files{i}=temp(i,1).name;
end
files=strrep(files,'.fa','');
files=strrep(files,'.hmm','');
files=strrep(files,'.out','');
files=strrep(files,'.faw','');
end
