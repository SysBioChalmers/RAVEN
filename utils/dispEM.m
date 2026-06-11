function dispEM(string,varargin)
% dispEM  Print a warning or throw an error.
%
% Helper function to print warnings or throw errors, optionally followed
% by a list of items.
%
% Parameters
% ----------
% string : char
%     the warning/error to show. "WARNING: " is appended automatically if
%     it is a warning.
% throwErrors : logical, optional
%     true if the function should throw an error (default true).
% toList : cell, optional
%     a cell array of items to list. If supplied, the string is printed
%     followed by each element in toList. If it is supplied but empty then
%     nothing is printed (default {}).
% trimWarnings : logical, optional
%     true if a maximum of 10 items should be displayed in a given
%     error/warning (default true).
%
% Examples
% --------
%     dispEM(string,throwErrors,toList,trimWarnings);

p=parseRAVENargs(varargin, {'throwErrors',true; 'toList','__notSupplied__'; 'trimWarnings',true});
throwErrors=p.throwErrors;
trimWarnings=p.trimWarnings;
toList=p.toList;
if ischar(toList) && isequal(toList,'__notSupplied__')
    toList=[];
elseif isempty(toList)
    return;
else
    toList=convertCharArray(toList);
end
if numel(toList)>10 && trimWarnings==true
    toList{10}=['...and ' num2str(numel(toList)-9) ' more'];
    toList(11:end)=[];
end
if throwErrors==false
    errorText=['WARNING: ' string '\n'];
    % Wrap text to command window size
    sz = get(0, 'CommandWindowSize');
    errorText = textwrap({errorText},sz(1));
    errorText = strjoin(errorText,'\n');
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
