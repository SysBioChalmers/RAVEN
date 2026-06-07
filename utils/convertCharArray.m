function inputConverted = convertCharArray(funcInput)
%convertCharArray
%   Converts input to make sure it is a cell array of character vectors.
%   String arrays are transformed into character vectors, and if only one
%   character vector is given. Output is always a cell array, also if only
%   one character vector is given as input.
%
%   Input:
%   funcInput          function input that should be checked
%
%   Output:
%   inputConverted      cell array of character vectors
%

if ~isempty(funcInput)
    if ischar(funcInput) || isstring(funcInput)
        inputConverted = {char(funcInput)};
    elseif ~iscell(funcInput) && ~(all(cellfun(@ischar, funcInput)) || all(cellfun(@isstring, funcInput)))
        error([inputname(1) ' should be a (cell array of) character vector(s)'])
    else
        inputConverted = cellstr(funcInput);
    end
else
    inputConverted=[];
end
end

