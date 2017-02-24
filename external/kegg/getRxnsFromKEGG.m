function model=getRxnsFromKEGG(keggPath,keepUndefinedStoich,keepIncomplete, keepGeneral)
% getRxnsFromKEGG
%   Loads keggRxns.mat
%
%   keggPath            path to the location where the KEGG files are stored,
%                       including keggRxns.mat
%   keepUndefinedStoich include reactions in the form n A <=> n+1 A. These
%                       will be dealt with as two separate metabolites
%                       (opt, default true)
%   keepIncomplete      include reactions which have been labelled as
%                       "incomplete", "erroneous" or "unclear" (opt,
%                       default true)
%   keepGeneral         include reactions which have been labelled as
%                       "general reaction". These are reactions on the form
%                       "an aldehyde <=> an alcohol", and are therefore
%                       unsuited for modelling purposes. Note that not all
%                       reactions have this type of annotation, and the
%                       script will therefore not be able to remove all
%                       such reactions (opt, default false)
%
%   model     a model structure generated from the database. The following
%             fields are filled
%             id:             'KEGG'
%             description:    'Automatically generated from KEGG database'
%             rxns:           KEGG reaction ids
%             rxnNames:       Name for each reaction entry
%             mets:           KEGG compound ids. If the equations use
%                             stoichiometry such as ID(n+1) then the whole
%                             expression is saved as the id
%             eccodes:        Corresponding ec-number if available
%             rxnMiriams:     Contains reaction specific information such as
%                             KO id and pathways that the reaction is
%                             associated to
%             S:              Stoichiometric matrix
%             lb:             -1000 for all reactions
%             ub:             1000 for all reactions
%             rev:            1 for reversible and 0 for irreversible. For
%                             reactions present in pathway maps the reversibility
%                             is taken from there
%             b:              0 for all metabolites
%
%   Reactions on the form A <=> A + B will not be loaded. This function no longer
%   supports generating new keggRxns.mat files from a KEGG FTP dump, as this
%   option has become obsolete.
%
%   Usage: model=getRxnsFromKEGG(keggPath,keepUndefinedStoich,keepIncomplete,keepGeneral)
%
%   Eduard Kerkhoven, 2017-02-24
%

if nargin<2
    keepUndefinedStoich=true;
end
if nargin<3
    keepIncomplete=true;
end
if nargin<4
    keepGeneral=false;
end

%Check if the reactions have been parsed before and saved. If so, load the
%model.
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
rxnsFile=fullfile(ravenPath,'external','kegg','keggRxns.mat');
fprintf(['NOTE: Importing KEGG reactions from ' strrep(rxnsFile,'\','/') '.\n']);
load(rxnsFile);

%Delete reaction which are labeled as "incomplete", "erroneous", "unclear"
%or "general reaction" (depending on settings.
if keepGeneral==false
    model=removeReactions(model,intersect(isGeneral,model.rxns),true,true);
end
if keepIncomplete==false
    model=removeReactions(model,intersect(isIncomplete,model.rxns),true,true);
end

%Delete reactions involving undefined stoichiometry. These metabolites have
%an ID containing the letter "n" or "m"
if keepUndefinedStoich==false
    I=cellfun(@any,strfind(model.mets,'n')) | cellfun(@any,strfind(model.mets,'m'));
    [~, J]=find(model.S(I,:));
    model=removeReactions(model,J,true,true);
end
end
