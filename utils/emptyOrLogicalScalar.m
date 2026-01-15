function emptyOrLogical(x)
    if isempty(x)
        return;
    end
    validateattributes(x, {'logical'}, {'scalar'});
end
