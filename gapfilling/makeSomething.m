function [solution, metabolite]=makeSomething(model,varargin)
% makeSomething  Excrete any metabolite using as few reactions as possible.
%
% Wrapper for findLeakMetabolite(model,'produce',...), retained for
% backward compatibility. See findLeakMetabolite for full documentation.
%
% See Also
% --------
% findLeakMetabolite
[solution, metabolite]=findLeakMetabolite(model,'produce',varargin{:});
end
