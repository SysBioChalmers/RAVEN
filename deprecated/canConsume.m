function consumed=canConsume(model,varargin)
% canConsume  DEPRECATED. Use canExchange(model,'consume',...).
%
% Deprecated wrapper kept so existing scripts keep working; it will be
% removed in the next major release.
%
% See Also
% --------
% canExchange

deprecationWarning('canConsume','canExchange(model,''consume'',...)');

consumed=canExchange(model,'consume',varargin{:});
end
