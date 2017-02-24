function model=getGenesFromKEGG(keggPath,koList)
% getGenesFromKEGG
%   Loads keggGenes.mat
%
%   keggPath    path to the location where the KEGG files are stored,
%               including keggGenes.mat
%
%   koList      the number of genes in KEGG is very large. koList can be a
%               cell array with KO identifiers, in which case only genes
%               belonging to one of those KEGG orthologies are retrieved (opt,
%               default all KOs with associated reactions)
%
%   model       a model structure generated from the database. The following
%               fields are filled
%               id:             'KEGG'
%               description:    'Automatically generated from KEGG database'
%               rxns:           KO ids
%               rxnNames:       Name for each entry
%               genes:          IDs for all the genes. Genes are saved as
%                               organism abbreviation:id (same as in KEGG).
%                               'HSA:124' for example is alcohol dehydrogenase
%                               in Homo sapiens
%               rxnGeneMat      A binary matrix that indicates whether a
%                               specific gene is present in a KO id
%
%   This function no longer supports generating new keggGenes.mat files from
%   a KEGG FTP dump, as this option has become obsolete.
%
%   Usage: model=getGenesFromKEGG(keggPath,koList)
%
%   Eduard Kerkhoven, 2017-02-24
%

[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
genesFile=fullfile(ravenPath,'external','kegg','keggGenes.mat');
fprintf(['NOTE: Importing KEGG genes from ' strrep(genesFile,'\','/') '.\n']);
load(genesFile);

%Only get the KOs in koList
I=~ismember(model.rxns,koList);
model=removeReactions(model,I,true,true);
end

function allKOs=getAllKOs(keggPath)
    %Retrieves all KOs that are associated to reactions. This is because
    %the number of genes in KEGG is very large so without this parsing it
    %would take many hours
    allKOs={};

    %First check if the reactions have already been parsed
    [ST I]=dbstack('-completenames');
    ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
    rxnsFile=fullfile(ravenPath,'external','kegg','keggRxns.mat');
    fprintf(['NOTE: Importing KEGG ORTHOLOGY list from ' strrep(rxnsFile,'\','/') '.\n']);
    load(rxnsFile,'model');
    %Loop through the reactions and add the corresponding genes
    for i=1:numel(model.rxns)
        if isstruct(model.rxnMiriams{i})
            %Get all KOs
            allKOs=[allKOs;model.rxnMiriams{i}.value(strcmpi(model.rxnMiriams{i}.name,'urn:miriam:kegg.ko'))];
        end
    end
    allKOs=unique(allKOs);
end
