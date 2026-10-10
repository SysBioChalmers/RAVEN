function out = parseYAML(filename)
% parseYAML  Read a YAML file into a MATLAB struct/cell tree.
%
% Supports block-style mappings and sequences (indentation-based),
% flow-style mappings and sequences (`{ key: value }`, `[a, b]`), quoted
% and unquoted scalars, and comments. Does not support multi-document
% files, anchors/aliases, tags, or block scalars (`|`, `>`).
%
% Use this for parsing YAML configuration / data files with a nested
% mapping/sequence/scalar shape (e.g. condition files consumed by
% applyCondition). For loading a cobra-format model YAML, use
% readYAMLmodel instead — that function knows the model schema and
% returns a populated RAVEN model struct.
%
% Parameters
% ----------
% filename : char
%     Path to the YAML file.
%
% Returns
% -------
% out : struct or cell or char or double or logical or []
%     MATLAB representation of the document:
%
%     - mapping  -> struct (keys sanitised to valid MATLAB field names)
%     - sequence -> cell column vector
%     - quoted or unquoted string -> char
%     - integer or float -> double
%     - true/false -> logical
%     - null / ~ / empty value -> []
%
% Examples
% --------
%     cfg = parseYAML('data/conditions/anaerobic.yml');
%
% See also
% --------
% readYAMLmodel

if ~isfile(filename)
    error('parseYAML:fileNotFound', 'File not found: %s', filename);
end

[lines, indents] = readYAMLLines(filename);
if isempty(lines)
    out = struct();
    return
end
[out, ~] = parseNode(lines, 1, indents);
end


function [lines, indents] = readYAMLLines(filename)
% Read filename into non-blank, comment-stripped content lines, paired
% with each line's leading-space indentation count. Document markers
% (`---`, `...`) are dropped.
fid = fopen(filename, 'r');
cleanup = onCleanup(@() fclose(fid));

lines = {};
indents = [];
tline = fgetl(fid);
while ischar(tline)
    content = stripComment(tline);
    trimmed = strtrim(content);
    if ~isempty(trimmed) && ~strcmp(trimmed, '---') && ~strcmp(trimmed, '...')
        indents(end+1,1) = numel(content) - numel(regexprep(content, '^ +', '')); %#ok<AGROW>
        lines{end+1,1} = trimmed; %#ok<AGROW>
    end
    tline = fgetl(fid);
end
end


function line = stripComment(line)
% Strip a trailing `# comment`, ignoring `#` inside a quoted string.
inSingle = false;
inDouble = false;
for i = 1:numel(line)
    c = line(i);
    if c == '''' && ~inDouble
        inSingle = ~inSingle;
    elseif c == '"' && ~inSingle
        inDouble = ~inDouble;
    elseif c == '#' && ~inSingle && ~inDouble
        line = line(1:i-1);
        break
    end
end
end


function [value, nextIdx] = parseNode(lines, idx, indents)
% Dispatch to a block sequence or block mapping, based on whether the
% line at idx starts a sequence item (`-`) or a mapping key.
if startsWithDash(lines{idx})
    [value, nextIdx] = parseSequence(lines, idx, indents);
else
    [value, nextIdx] = parseMapping(lines, idx, indents);
end
end


function tf = startsWithDash(s)
tf = ~isempty(s) && (strcmp(s, '-') || (numel(s) >= 2 && s(1) == '-' && s(2) == ' '));
end


function [value, nextIdx] = parseMapping(lines, idx, indents)
baseIndent = indents(idx);
value = struct();
while idx <= numel(lines) && indents(idx) == baseIndent && ~startsWithDash(lines{idx})
    [key, rest] = splitKeyValue(lines{idx});
    key = matlabFieldName(key);
    if isempty(rest)
        if idx+1 <= numel(lines) && indents(idx+1) > baseIndent
            [childVal, idx] = parseNode(lines, idx+1, indents);
        else
            childVal = [];
            idx = idx + 1;
        end
    else
        childVal = parseScalarOrFlow(rest);
        idx = idx + 1;
    end
    value.(key) = childVal;
end
nextIdx = idx;
end


function [value, nextIdx] = parseSequence(lines, idx, indents)
baseIndent = indents(idx);
items = {};
while idx <= numel(lines) && indents(idx) == baseIndent && startsWithDash(lines{idx})
    rest = strtrim(lines{idx}(2:end));
    if isempty(rest)
        if idx+1 <= numel(lines) && indents(idx+1) > baseIndent
            [itemVal, idx] = parseNode(lines, idx+1, indents);
        else
            itemVal = [];
            idx = idx + 1;
        end
    else
        itemVal = parseScalarOrFlow(rest);
        idx = idx + 1;
    end
    items{end+1,1} = itemVal; %#ok<AGROW>
end
value = items;
nextIdx = idx;
end


function [key, rest] = splitKeyValue(line)
% Split "key: value" (or "key:") on the first top-level colon.
inSingle = false;
inDouble = false;
colonPos = 0;
for i = 1:numel(line)
    c = line(i);
    if c == '''' && ~inDouble
        inSingle = ~inSingle;
    elseif c == '"' && ~inSingle
        inDouble = ~inDouble;
    elseif c == ':' && ~inSingle && ~inDouble
        colonPos = i;
        break
    end
end
if colonPos == 0
    error('parseYAML:invalidLine', 'Expected "key: value" or "key:", got: %s', line);
end
key = strtrim(line(1:colonPos-1));
key = strrep(strrep(key, '"', ''), '''', '');
rest = strtrim(line(colonPos+1:end));
end


function value = parseScalarOrFlow(s)
s = strtrim(s);
if ~isempty(s) && s(1) == '{'
    value = parseFlowMap(s);
elseif ~isempty(s) && s(1) == '['
    value = parseFlowSeq(s);
else
    value = parseScalar(s);
end
end


function value = parseFlowMap(s)
inner = strtrim(s(2:end-1));
value = struct();
if isempty(inner)
    return
end
parts = splitTopLevel(inner, ',');
for i = 1:numel(parts)
    [key, rest] = splitKeyValue(parts{i});
    value.(matlabFieldName(key)) = parseScalarOrFlow(rest);
end
end


function value = parseFlowSeq(s)
inner = strtrim(s(2:end-1));
if isempty(inner)
    value = cell(0,1);
    return
end
parts = splitTopLevel(inner, ',');
value = cell(numel(parts), 1);
for i = 1:numel(parts)
    value{i} = parseScalarOrFlow(parts{i});
end
end


function parts = splitTopLevel(s, delim)
% Split s on delim, ignoring delim inside quotes or nested {}/[].
parts = {};
depth = 0;
inSingle = false;
inDouble = false;
last = 1;
for i = 1:numel(s)
    c = s(i);
    if c == '''' && ~inDouble
        inSingle = ~inSingle;
    elseif c == '"' && ~inSingle
        inDouble = ~inDouble;
    elseif ~inSingle && ~inDouble
        if c == '{' || c == '['
            depth = depth + 1;
        elseif c == '}' || c == ']'
            depth = depth - 1;
        elseif c == delim && depth == 0
            parts{end+1} = s(last:i-1); %#ok<AGROW>
            last = i + 1;
        end
    end
end
parts{end+1} = s(last:end); %#ok<AGROW>
parts = strtrim(parts);
end


function value = parseScalar(s)
if isempty(s) || strcmp(s, '~') || strcmpi(s, 'null')
    value = [];
    return
end
if numel(s) >= 2 && ((s(1) == '"' && s(end) == '"') || (s(1) == '''' && s(end) == ''''))
    value = s(2:end-1);
    return
end
if strcmpi(s, 'true')
    value = true;
    return
end
if strcmpi(s, 'false')
    value = false;
    return
end
if ~isempty(regexp(s, '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$', 'once'))
    value = str2double(s);
    return
end
value = s;
end


function name = matlabFieldName(key)
% Sanitise a YAML key into a valid MATLAB field name. Replaces non-
% alphanumeric characters with underscores; prefixes a digit-starting
% key with 'f_'.
name = regexprep(key, '[^A-Za-z0-9_]', '_');
if isempty(name) || ~isstrprop(name(1), 'alpha')
    name = ['f_' name];
end
end
