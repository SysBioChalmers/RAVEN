function arrayData=parseHPArna(fileName, version)
% parseHPA
%   Parses a database dump of the Human Protein Atlas (HPA) RNA-Seq data.
%
%   Input:
%   fileName            tab-separated database dump of HPA RNA data. For
%                       details regarding the format, see
%                       http://www.proteinatlas.org/about/download.
%   version             version of HPA [optional, default=18]
%
%
%   Output:
%   arrayData
%       genes               cell array with the unique ensemble gene IDs
%       geneNames           cell array with the gene names (gene abbrevs)
%       tissues             cell array with the tissue names
%       levels              matrix of gene expression levels (TPM), where
%                           rows correspond to genes, and columns
%                           correspond to tissues
%
%   Usage: arrayData=parseHPArna(fileName,version)

if nargin<2
    %Change this and add code for more versions when the current HPA
    %version is increased and the format is changed
    version=18;
end

if ~(exist(fileName,'file')==2)
    error('HPA file %s cannot be found',string(fileName));
end

if (version == 18)
    fid=fopen(fileName,'r');
    hpa=textscan(fid,'%q %q %q %q %q','Delimiter','\t');
    fclose(fid);
    
    %Go through and see if the headers match what was expected
    headers={'Gene' 'Gene name' 'Sample' 'Value' 'Unit'};
    for i=1:numel(headers)
        if ~strcmpi(headers(i),hpa{i}(1))
            EM=['Could not find the header "' headers{i} '". Make sure that the input file matches the format specified at http://www.proteinatlas.org/about/download'];
            dispEM(EM);
        end
        %Remove the header line here
        hpa{i}(1)=[];
    end
    
    %Get unique gene IDs and tissue names
    [arrayData.genes, P, I] = unique(hpa{1});
    arrayData.geneNames = hpa{2}(P);  % retrieve corresponding gene names
    [arrayData.tissues, ~, J] = unique(hpa{3});
    
    %Now extract the gene levels and organize into matrix
    arrayData.levels = NaN(max(I),max(J));
    linearInd = sub2ind(size(arrayData.levels),I,J);
    arrayData.levels(linearInd) = str2double(hpa{4});
    
else
    error('Previous HPA versions (before 18) not currently supported.');
end
