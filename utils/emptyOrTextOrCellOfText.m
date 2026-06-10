function mustBeEmptyOrTextOrCellOfText(x)
% mustBeEmptyOrTextOrCellOfText  Validate empty, text scalar or cell of text.
%
% Validate [] OR a text scalar (char row or string scalar) OR a cell
% array of such text. Throws an error if the input is none of these.
%
% Parameters
% ----------
% x : double or char or string or cell
%     value to validate; allowed forms are empty ([]), a character row
%     vector, a string scalar, or a cell array whose elements are each a
%     character row vector or string scalar.

% Allow explicit empty
if isempty(x)
    return;
end

% Allow a single text scalar
if ischar(x)
    if ~isrow(x)
        error('Text must be a character row vector, a string scalar, or empty.');
    end
    return;
elseif isstring(x)
    if ~isscalar(x)
        error('Text must be a string scalar, a character row vector, or empty.');
    end
    return;
end

% Allow a cell array of text
if iscell(x)
    for k = 1:numel(x)
        v = x{k};
        if ischar(v)
            if ~isrow(v)   % '' is 1x0 char, which is a row -> allowed
                error('Each cell element must be a character row vector or a string scalar.');
            end
        elseif isstring(v)
            if ~isscalar(v)  % "" is a scalar string -> allowed
                error('Each cell element must be a string scalar or a character row vector.');
            end
        else
            error('Each cell element must be text (string scalar or character row vector).');
        end
    end
    return;
end

% Anything else is invalid
error('Value must be empty, a text value, or a cell array of text values.');
end
