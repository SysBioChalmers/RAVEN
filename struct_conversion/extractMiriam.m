function [miriams,extractedMiriamNames]=extractMiriam(modelMiriams,miriamNames)
% extractMiriam
%   This function unpacks the information kept in metMiriams, rxnMiriams,
%   geneMiriams or compMiriams to make the annotation more
%   human-readable. The obtained cell array looks the same like in Excel
%   format, just the columns are split to have particular miriam name in
%   corresponding column
%
%   modelMiriams                a miriam structure (e.g. model.metMiriams)
%                               for one or multiple metabolites
%   miriamNames                 cell array with miriam names to be
%                               extracted (optional, default 'all', meaning
%                               that annotation for all miriam names found
%                               in modelMiriams will be extracted)
%
%   miriams                     a cell array with extracted miriams. if
%                               several miriam names are requested, the
%                               corresponding information is saved in
%                               different columns. if there are several ids
%                               available for the same entity (metabolite,
%                               gene, reaction or compartment), they are
%                               concatenated into one column. the total
%                               number of column represent the number of
%                               unique miriam names per entity
%   extractedMiriamNames        cell array with extracted miriam names
%
%   Usage: miriam=extractMiriam(modelMiriams,miriamName,addNull)
%
%   Simonas Marcisauskas, 2018-04-12
%   Benjamin Sanchez, 2018-04-13
%

if nargin<2 || strcmp(miriamNames,'all')
    extractAllTypes=true;
else
    extractAllTypes=false;
end

%The annotation for all miriam names should be extracted
if extractAllTypes
    miriamNames='';
    for i=1:numel(modelMiriams)
        if ~isempty(modelMiriams{i,1})
            for j=1:numel(modelMiriams{i,1}.name)
                miriamNames{numel(miriamNames)+1,1}=modelMiriams{i,1}.name{j};
            end
        end
    end
    miriamNames=sort(unique(miriamNames));
end

%Ensure that the list of miriam names is in appropriate structure
if ischar(miriamNames)
    miriamNames={miriamNames};
end

%Aggregate the final cell array table with extracted miriam information
miriams='';
for i=1:numel(miriamNames)
    miriams=[miriams, extractMiriamType(modelMiriams,miriamNames{i})];
end
extractedMiriamNames=miriamNames;

end

function miriams=extractMiriamType(modelMiriams,miriamName)
%Create an empty cell array for ids vector
tempMiriams = cell([size(modelMiriams,1) 1]);
%Firstly obtain the list of relevant miriam ids. Several entries may have
%several miriam ids, such ids are kept in additional columns
for i=1:numel(modelMiriams)
    if (~isempty(modelMiriams{i,1})) && any(strcmp(modelMiriams{i,1}.name,miriamName))
        for j=1:numel(modelMiriams{i,1}.name)
            if strcmp(modelMiriams{i,1}.name(j),miriamName)
                tempMiriams(i,j) = modelMiriams{i,1}.value(j);
            end
        end
    else
        modelMiriams{i,1} = '';
    end
end

%Now add miriam names to the newly obtained miriam ids. Entities with
%multiple ids per miriam name have the ids concatenated by using semicolon
%as separator
miriams = cell([size(tempMiriams,1) 1]);
for i=1:size(tempMiriams,1)
    for j=1:size(tempMiriams,2)
        if j==1
            miriams{i,1}=strcat(miriamName,'/',tempMiriams{i,1});
        else
            miriams{i,1}=strcat(miriams{i,1},';',miriamName,'/',tempMiriams{i,j});
        end
    end
end

%Ensure that cell positions without miriams are blank
miriams=regexprep(miriams,strcat(miriamName,'/;'),'');
miriams=regexprep(miriams,strcat('^',miriamName,'/$'),'');
miriams=regexprep(miriams,strcat(';',miriamName,'/$'),'');

%Delete middle names:
miriams=regexprep(miriams,strcat(';',miriamName,'/'),';');

end
