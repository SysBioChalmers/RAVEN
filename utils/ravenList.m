function msg = ravenList(header, items, trim)
% ravenList  Format a header and item list for warning() or error().
%
% Parameters
% ----------
% header : char
%     message header string.
% items : cell
%     cell array of items to list below the header.
% trim : logical
%     cap the list at 10 items (default true).
%
% Returns
% -------
% msg : char
%     header + newline + tab-indented items, ready to pass to warning() or
%     error() as the format string.
%
% Examples
% --------
%     warning('RAVEN:warning', '%s', ravenList('Missing metabolites:', badMets));
%     error('RAVEN:badInput', '%s', ravenList('Invalid reactions:', badRxns, false));

if nargin < 3; trim = true; end
items = convertCharArray(items);
if trim && numel(items) > 10
    items{10} = sprintf('...and %d more', numel(items)-9);
    items(11:end) = [];
end
msg = [header, newline, strjoin(strcat(char(9), items), newline)];
end
