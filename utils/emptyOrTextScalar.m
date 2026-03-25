function mustBeEmptyOrTextScalar(x)
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
