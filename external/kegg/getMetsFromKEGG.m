function model=getMetsFromKEGG(keggPath)
% getMetsFromKEGG
%   Loads keggMets.mat
%
%   keggPath    path to the location where the KEGG files are stored,
%               including keggMets.mat
%
%   model       a model structure generated from the database. The following
%               fields are filled
%               id:             'KEGG'
%               description:    'Automatically generated from KEGG database'
%               mets:           KEGG compound ids
%               metNames:       Compound name. Only the first name will be
%                               saved if there are several synonyms
%               metMiriams:     If there is a CHEBI id available, then that
%                               will be saved here
%               inchis:         InChI string for the metabolite
%               metFormulas:    The chemical composition of the metabolite.
%                               This will only be loaded if there is no InChI 
%                               string
%
%   This function no longer supports generating new keggMets.mat files from
%   a KEGG FTP dump, as this option has become obsolete.
%
%   Usage: model=getMetsFromKEGG(keggPath)
%
%   Eduard Kerkhoven, 2017-02-24
%

[ST I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
metsFile=fullfile(ravenPath,'external','kegg','keggMets.mat');
fprintf(['NOTE: Importing KEGG metabolites from ' strrep(metsFile,'\','/') '.\n']);
load(metsFile);
