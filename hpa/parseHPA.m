function hpaData=parseHPA(fileName, version)
% parseHPA
%   Parses a database dump of the Human Protein Atlas (HPA)
%
%   fileName            comma- or tab-separated database dump of HPA. For details
%                       regarding the format, see
%                       http://www.proteinatlas.org/about/download.
%   version             version of HPA [optional, default=17]
%
%   hpaData
%       genes               cell array with the unique gene names. In
%                           version 17 this is the ensamble name, see
%                           geneNames below for the names in ver 17
%       geneNames           cell array with the gene names, indexed the
%                           same way as genes.
%       tissues             cell array with the tissue names. The list may not be
%                           unique, as there can be multiple cell types per tissue
%       celltypes           cell array with the cell type names for each tissue
%       levels              cell array with the unique expression levels
%       types               cell array with the unique evidence types
%       reliabilities       cell array with the unique reliability levels
%
%       gene2Level          gene-to-expression level mapping in sparse matrix form.
%                           The value for element i,j is the index in
%                           hpaData.levels of gene i in cell type j
%       gene2Type           gene-to-evidence type mapping in sparse matrix form.
%                           The value for element i,j is the index in
%                           hpaData.types of gene i in cell type j. Doesn't
%                           exist in version 17.
%       gene2Reliability    gene-to-reliability level mapping in sparse matrix form.
%                           The value for element i,j is the index in
%                           hpaData.reliabilities of gene i in cell type j
%
%
%   Usage: hpaData=parseHPA(fileName)
%
%   Johan Gustafsson, 2017-10-10

if nargin<2
    version=17; %Change this and add code for more versions when the current HPA version is increased and the format is changed
end;

if (version == 17)
    fid=fopen(fileName,'r');
    hpa=textscan(fid,'%q %q %q %q %q %q','Delimiter','\t');
    fclose(fid);

    %Go through and see if the headers match what was expected
    headers={'Gene' 'Gene name' 'Tissue' 'Cell type' 'Level' 'Reliability'};
    for i=1:numel(headers)
        if ~strcmpi(headers(i),hpa{i}(1))
            EM=['Could not find the header "' headers{i} '". Make sure that the input file matches the format specified at http://www.proteinatlas.org/about/download'];
          dispEM(EM);
        end
        %Remove the header line here
        hpa{i}(1)=[];
    end

    %Get the unique values of each data type
    [hpaData.genes, P, I]=unique(hpa{1});
    hpaData.geneNames=hpa{2}(P); %make this vector use the index as genes
    [~, J, K]=unique(strcat(hpa{3},'€',hpa{4}));
    hpaData.tissues=hpa{3}(J);
    hpaData.celltypes=hpa{4}(J);
    [hpaData.levels, ~, L]=unique(hpa{5});
    [hpaData.reliabilities, ~, N]=unique(hpa{6});

    %Map the data to be sparse matrises instead
    hpaData.gene2Level=sparse(I,K,L,numel(hpaData.genes),numel(hpaData.tissues));
    hpaData.gene2Reliability=sparse(I,K,N,numel(hpaData.genes),numel(hpaData.tissues));
else
    fid=fopen(fileName,'r');
    hpa=textscan(fid,'%q %q %q %q %q %q','Delimiter',',');
    fclose(fid);

    %Go through and see if the headers match what was expected
    headers={'Gene' 'Tissue' 'Cell type' 'Level' 'Expression type' 'Reliability'};
    for i=1:numel(headers)
        if ~strcmpi(headers(i),hpa{i}(1))
            EM=['Could not find the header "' headers{i} '". Make sure that the input file matches the format specified at http://www.proteinatlas.org/about/download'];
          dispEM(EM);
        end
        %Remove the header line here
        hpa{i}(1)=[];
    end

    %Get the unique values of each data type
    [hpaData.genes, ~, I]=unique(hpa{1});
    [~, J, K]=unique(strcat(hpa{2},'€',hpa{3}));
    hpaData.tissues=hpa{2}(J);
    hpaData.celltypes=hpa{3}(J);
    [hpaData.levels, ~, L]=unique(hpa{4});
    [hpaData.types, ~, M]=unique(hpa{5});
    [hpaData.reliabilities, ~, N]=unique(hpa{6});

    %Map the data to be sparse matrises instead
    hpaData.gene2Level=sparse(I,K,L,numel(hpaData.genes),numel(hpaData.tissues));
    hpaData.gene2Type=sparse(I,K,M,numel(hpaData.genes),numel(hpaData.tissues));
    hpaData.gene2Reliability=sparse(I,K,N,numel(hpaData.genes),numel(hpaData.tissues));    
end