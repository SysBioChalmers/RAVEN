function mustBeEmptyOrTextScalar(x)
% emptyOrTextScalar  Validate that a value is empty or a string/char scalar.
%
% Argument-validation helper used in arguments blocks. Passes silently if x
% is empty, a string scalar, or a char row vector, and otherwise raises a
% validation error.
%
% Parameters
% ----------
% x : char or string
%     the value to validate; may be empty.
    if isempty(x)
        return;
    end
    if isstring(x)
        if ~isscalar(x), error('Value must be a string scalar or empty.'); end
    elseif ischar(x)
        if ~isrow(x), error('Char input must be a row vector or empty.'); end
    else
        error('Value must be a string/char scalar or empty.');
    end
end
