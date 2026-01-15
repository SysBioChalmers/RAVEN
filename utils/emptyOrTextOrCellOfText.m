function mustBeEmptyOrTextOrCellOfText(x)
% Validate [] OR text scalar (char row or string scalar) OR cell array of such text.

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
