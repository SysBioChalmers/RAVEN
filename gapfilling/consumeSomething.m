function [solution, metabolite]=consumeSomething(model,varargin)
% consumeSomething  Consume any metabolite using as few reactions as possible.
%
% Wrapper for findLeakMetabolite(model,'consume',...), retained for
% backward compatibility. See findLeakMetabolite for full documentation.
%
% See Also
% --------
% findLeakMetabolite
[solution, metabolite]=findLeakMetabolite(model,'consume',varargin{:});
end
