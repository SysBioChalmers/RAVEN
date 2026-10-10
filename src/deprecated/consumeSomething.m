function [solution, metabolite]=consumeSomething(model,varargin)
% consumeSomething  DEPRECATED. Use findLeakMetabolite(model,'consume',...).
%
% Deprecated wrapper kept so existing scripts keep working; it will be
% removed in the next major release.
%
% Its positional argument order predates findLeakMetabolite and is
% preserved here:
%
%     consumeSomething(model, ignoreMets, isNames, minNrFluxes, params, ignoreIntBounds)
%
% findLeakMetabolite inserts allowExcretion (a produce-only option) in
% position 4, so the arguments are translated rather than forwarded
% verbatim. Forwarding them verbatim, as this wrapper previously did,
% silently shifted params into allowExcretion and ignoreIntBounds into
% params, so a requested ignoreIntBounds was dropped.
%
% See Also
% --------
% findLeakMetabolite

deprecationWarning('consumeSomething','findLeakMetabolite(model,''consume'',...)');

p=parseRAVENargs(varargin, {'ignoreMets',[]; ...
    'isNames',false; ...
    'minNrFluxes',false; ...
    'params',[]; ...
    'ignoreIntBounds',false});

[solution, metabolite]=findLeakMetabolite(model,'consume', ...
    'ignoreMets',p.ignoreMets, ...
    'isNames',p.isNames, ...
    'minNrFluxes',p.minNrFluxes, ...
    'params',p.params, ...
    'ignoreIntBounds',p.ignoreIntBounds);
end
