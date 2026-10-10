function [solution, metabolite]=makeSomething(model,varargin)
% makeSomething  DEPRECATED. Use findLeakMetabolite(model,'produce',...).
%
% Deprecated wrapper kept so existing scripts keep working; it will be
% removed in the next major release. Its argument order already matches
% findLeakMetabolite's, so the arguments are forwarded as they are.
%
% See Also
% --------
% findLeakMetabolite

deprecationWarning('makeSomething','findLeakMetabolite(model,''produce'',...)');

[solution, metabolite]=findLeakMetabolite(model,'produce',varargin{:});
end
