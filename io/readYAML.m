function out = readYAML(filename)
% readYAML
%   Read an arbitrary YAML file into a MATLAB struct / cell tree.
%
%   Use this for parsing arbitrary YAML configuration / data files
%   (e.g. yeast-GEM's data/conditions/*.yml). For loading a cobra-format
%   model YAML, use readYAMLmodel instead — that function knows the
%   model schema and returns a populated RAVEN model struct.
%
%   Implementation: delegates to Python's yaml.safe_load, then
%   recursively converts the py.dict / py.list tree to native MATLAB
%   struct / cell. Requires a working MATLAB-Python bridge and the
%   pyyaml package in the linked Python environment:
%
%       pip install pyyaml         % from the MATLAB-linked Python env
%
%   Input:
%       filename    path to the YAML file.
%
%   Output:
%       out         MATLAB representation of the document:
%                       py.dict   -> struct
%                       py.list   -> cell column vector
%                       py.str    -> char
%                       py.int    -> double
%                       py.float  -> double
%                       py.bool   -> logical
%                       py.None   -> []
%
% Usage: cfg = readYAML('data/conditions/anaerobic.yml')

if ~isfile(filename)
    error('readYAML:fileNotFound', 'File not found: %s', filename);
end

try
    py.importlib.import_module('yaml');
catch ME
    error('readYAML:pyyamlMissing', ...
        ['pyyaml is required to read arbitrary YAML files. Install it ' ...
         'in your MATLAB-linked Python environment (`pip install pyyaml`).' ...
         '\nUnderlying error: %s'], ME.message);
end

f = py.builtins.open(filename, 'r');
cleanup = onCleanup(@() f.close());
data = py.yaml.safe_load(f);

out = pyToMatlab(data);
end


function v = pyToMatlab(obj)
% Recursively convert pyyaml-loaded Python objects into MATLAB types.
if isa(obj, 'py.NoneType')
    v = [];
elseif isa(obj, 'py.bool')
    v = logical(obj);
elseif isa(obj, 'py.int') || isa(obj, 'py.float')
    v = double(obj);
elseif isa(obj, 'py.str')
    v = char(obj);
elseif isa(obj, 'py.dict')
    v = struct();
    keys = cell(py.list(obj.keys()));
    vals = cell(py.list(obj.values()));
    for i = 1:numel(keys)
        v.(matlabFieldName(char(keys{i}))) = pyToMatlab(vals{i});
    end
elseif isa(obj, 'py.list') || isa(obj, 'py.tuple')
    cells = cell(obj);
    v = cell(numel(cells), 1);
    for i = 1:numel(cells)
        v{i} = pyToMatlab(cells{i});
    end
else
    % Fallback: best-effort
    v = obj;
end
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
