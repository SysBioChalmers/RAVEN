function dispEM(string,throwErrors,toList,trimWarnings)
% dispEM
%   Helper function to print warning/errors
%
%   string          the warning/error to show. "WARNING: " is appended automatically
%                   if a warning
%   throwErrors     true if the function should throw an error (optional, default true)
%   toList          a cell array of items to list. If supplied, then the
%                   string will be printed followed by each element in
%                   toList. If it is supplied but empty then nothing is
%                   printed (optional, default {})
%   trimWarnings    true if only a maximal of 10 items should be displayed in
%                   a given error/warning (optional, default true)
%
% Usage: dispEM(string,throwErrors,toList,trimWarnings)

if nargin<2
    throwErrors=true;
end
if nargin<3
    toList=[];
elseif isempty(toList)
    return;
else
    toList=convertCharArray(toList);
end
if nargin<4
    trimWarnings=true;
end
if numel(toList)>10 && trimWarnings==true
    toList{10}=['...and ' num2str(numel(toList)-9) ' more'];
    toList(11:end)=[];
end
if throwErrors==false
    errorText=['WARNING: ' string '\n'];
else
    errorText=[string '\n'];
end
if ~isempty(toList)
    for i=1:numel(toList)
        errorText=[errorText '\t' toList{i} '\n'];
    end
end
if throwErrors==false
    %Escape special characters, required for fprintf
    errorText=regexprep(errorText,'(\\|\%|'')(\\n)$','\\$0');
    fprintf([errorText '\n']);
else
    throw(MException('',errorText));
end
end
