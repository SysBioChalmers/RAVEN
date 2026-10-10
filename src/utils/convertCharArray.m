function inputConverted = convertCharArray(funcInput)
% convertCharArray  Ensure input is a cell array of character vectors.
%
% Converts input to make sure it is a cell array of character vectors.
% String arrays are transformed into character vectors. Output is always a
% cell array, also if only one character vector is given as input.
%
% Parameters
% ----------
% funcInput : char or cell or string
%     function input that should be checked.
%
% Returns
% -------
% inputConverted : cell
%     cell array of character vectors.
%
% Examples
% --------
%     inputConverted = convertCharArray(funcInput);

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

