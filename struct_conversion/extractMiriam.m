function miriams=extractMiriam(modelMiriams,miriamName,addNull)
% extractMiriam
%   This function is useful if miriams are needed to be extracted to the
%   separate cell array, e.g. ChEBI ids or KEGG COMPOUND ids
%
%   modelMiriams      a miriam structure (e.g. model.metMiriams)
%   miriamName        the name of miriam which is supposed to be extracted
%                     (e.g. ChEBI ids)
%   addNull           true if the indexes without miriams should contain
%                     'null' strings, otherwise corresponding positions are
%                     left blank (optional, default false)
%
%   miriams            a cell array with extracted miriams
%
%   Usage: miriam=extractMiriam(modelMiriams,miriamName,addNull)
%
%   Simonas Marcisauskas, 2017-06-01
%

if nargin<3
	addNull=false;
end;
if nargin<4
	header='';
end;

% Creating an empty cell array for rxn ids vector;
tempMiriams = cell([size(modelMiriams,1) 1]);
% Firstly obtaining the list of relevant miriam ids. Several entries
% may have several miriam ids, such ids are kept in additional columns;
for i=1:numel(modelMiriams)
    if (~isempty(modelMiriams{i,1})) && any(strcmp(modelMiriams{i,1}.name,miriamName))
        for j=1:numel(modelMiriams{i,1}.name)
            if strcmp(modelMiriams{i,1}.name(j),miriamName)
                tempMiriams(i,j) = modelMiriams{i,1}.value(j);
            end;           
        end;
    else
        modelMiriams{i,1} = 'null';
        if addNull==false
            tempMiriams{i,1} = regexprep(modelMiriams{i,1},'null','');
        end;
    end;
end;
% Now adding headers to the newly obtained miriam ids. For reactions
% which are associated with more than one miriam, these are
% concatenated by using semicolon as separator;
miriams = cell([size(tempMiriams,1) 1]);
for i=1:size(tempMiriams,1)
    for j=1:size(tempMiriams,2)
        if j==1
            miriams{i,1}=strcat(miriamName,'/',tempMiriams{i,1});
        else
            miriams{i,1}=strcat(miriams{i,1},';',miriamName,'/',tempMiriams{i,j});
        end;
    end;
    % Make sure that cell positions without miriams are blank;
    miriams{i,1}=regexprep(miriams{i,1},strcat(miriamName,'/;'),'');
    miriams{i,1}=regexprep(miriams{i,1},strcat('^',miriamName,'/$'),'');
end;

end