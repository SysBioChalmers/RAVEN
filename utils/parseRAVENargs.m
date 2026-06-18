function args = parseRAVENargs(rawArgs, spec)
% parseRAVENargs  Resolve optional arguments passed positionally or by name.
%
% Helper that lets a RAVEN function accept its optional arguments either in
% the traditional positional order, as name-value pairs, or as a mix of
% leading positional values followed by name-value pairs. The required
% leading arguments stay explicit in the function signature; the optional
% tail is collected into varargin and passed here together with a
% specification of the optional parameters.
%
% Parameters
% ----------
% rawArgs : cell
%     the varargin of the function, i.e. whatever was supplied after the required
%     positional arguments.
% spec : cell
%     an N-by-2 or N-by-3 cell array describing the optional parameters, one
%     per row:
%
%     - column 1 : parameter name (char).
%     - column 2 : default value, used when the parameter is not supplied.
%     - column 3 : optional validator, a function handle that is called with
%       the resolved value and should error on an invalid value (e.g. one of
%       the utils/emptyOr* helpers). Use [] for no validation.
%
% Returns
% -------
% args : struct
%     a struct with one field per parameter name, holding the resolved value.
%
% Notes
% -----
% Calling convention is resolved as follows:
%
%   1. If rawArgs is empty, all parameters take their defaults.
%   2. The first element in rawArgs that is a string matching a known
%      parameter name marks position k. Everything before k is treated
%      as positional (assigned left-to-right by the order in spec).
%      Everything from k onward must form valid name-value pairs: an even
%      count where every name-position element is a known parameter name.
%   3. If no such position exists (no element is a matching string), or the
%      trailing section from k does not form valid name-value pairs, the
%      entire rawArgs is treated as positional.
%
% This means three styles are accepted interchangeably:
%
%     f(model, v1, v2)                   % purely positional
%     f(model, "p1", v1, "p2", v2)       % purely named
%     f(model, v1, "p2", v2, "p3", v3)   % hybrid: v1 positional, rest named
%
% When a positional value could itself be a string that equals a parameter
% name, use the explicit name-value form for that argument to avoid
% ambiguity.
%
% Examples
% --------
%     % inside exportToExcelFormat(model, varargin):
%     p = parseRAVENargs(varargin, { ...
%         "fileName", []; ...
%         "sortIds",  false});
%     fileName = p.fileName;
%     sortIds  = p.sortIds;
%
%     % the caller may use any of these:
%     exportToExcelFormat(model, "model.xlsx", true)              % positional
%     exportToExcelFormat(model, "sortIds", true)                 % named
%     exportToExcelFormat(model, "model.xlsx", "sortIds", true)   % hybrid

if nargin < 2 || isempty(spec)
    error('parseRAVENargs:noSpec', 'A parameter specification is required.');
end

names = spec(:, 1);
defaults = spec(:, 2);
hasValidators = size(spec, 2) >= 3;

% Seed every parameter with its default.
args = struct();
for i = 1:numel(names)
    args.(names{i}) = defaults{i};
end

if ~isempty(rawArgs)
    % Find which elements are strings matching a known parameter name.
    isParamName = cellfun(@(x) (ischar(x) || isstring(x)) && isscalar(string(x)) && ...
                          any(strcmpi(x, names)), rawArgs);

    % Position of the first parameter-name string — where name-value begins.
    nvStart = find(isParamName, 1);

    nvOK = false;
    if ~isempty(nvStart)
        nv = rawArgs(nvStart:end);
        % Valid name-value section: even count, all name-position elements
        % are known parameter names.
        nvOK = mod(numel(nv), 2) == 0 && ...
               all(cellfun(@(x) (ischar(x) || isstring(x)) && isscalar(string(x)) && ...
                   any(strcmpi(x, names)), nv(1:2:end)));
    end

    if ~isempty(nvStart) && nvOK
        % Assign leading positional arguments (may be zero when nvStart == 1).
        nPos = nvStart - 1;
        if nPos > numel(names)
            error('parseRAVENargs:tooManyArgs', ...
                'Too many positional arguments: expected at most %d, got %d.', ...
                numel(names), nPos);
        end
        for k = 1:nPos
            args.(names{k}) = rawArgs{k};
        end
        % Assign name-value pairs.
        nv = rawArgs(nvStart:end);
        for k = 1:2:numel(nv)
            args.(names{strcmpi(nv{k}, names)}) = nv{k + 1};
        end
    else
        % Purely positional.
        if numel(rawArgs) > numel(names)
            error('parseRAVENargs:tooManyArgs', ...
                'Too many positional arguments: expected at most %d, got %d.', ...
                numel(names), numel(rawArgs));
        end
        for k = 1:numel(rawArgs)
            args.(names{k}) = rawArgs{k};
        end
    end
end

% Apply optional validators.
if hasValidators
    for i = 1:numel(names)
        validator = spec{i, 3};
        if ~isempty(validator)
            validator(args.(names{i}));
        end
    end
end
end
