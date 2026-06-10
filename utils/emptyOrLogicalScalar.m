function emptyOrLogical(x)
% emptyOrLogicalScalar  Validate that a value is empty or a logical scalar.
%
% Argument-validation helper used in arguments blocks. Passes silently if x
% is empty or a logical scalar, and otherwise raises a validation error.
%
% Parameters
% ----------
% x : logical
%     the value to validate; may be empty.
    if isempty(x)
        return;
    end
    validateattributes(x, {'logical'}, {'scalar'});
end
