function phylDistStruct=getPhylDist(keggPath,onlyInKingdom)
% getPhylDist
%   Loads keggPhylDist.mat
%
%   keggPath        path to the location where the KEGG files are stored,
%                   including keggPhylDist.mat
%   onlyInKingdom   if true, it generates a distance matrix with distance 
%                   Inf for organisms from another domains
%                   (Prokaryota, Eukaryota) (opt, default false)
%
%   phylDistStruct  a structure with a list of organism ids and a matrix
%                   that specifies their pairwise distances
%
%   This simple metric is based on the number of nodes two organisms are
%   away from each other in KEGG. This function no longer supports
%   generating new keggPhylDist.mat files from a KEGG FTP dump, as this
%   option has become obsolete.
%
%   Usage: phylDistStruct=getPhylDist(keggPath,onlyInKingdom)
%
%   Eduard Kerkhoven, 2017-02-24
%

if nargin<2
    onlyInKingdom=false;
end

%Check if the reactions have been parsed before and saved. If so, load the
%model.
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
distFile=fullfile(ravenPath,'external','kegg','keggPhylDist.mat');
fprintf(['NOTE: Importing KEGG phylogenetic distance matrix from ' strrep(distFile,'\','/') '.\n']);
load(distFile);
