function produced=canProduce(model,varargin)
% canProduce  DEPRECATED. Use canExchange(model,'produce',...).
%
% Deprecated wrapper kept so existing scripts keep working; it will be
% removed in the next major release.
%
% See Also
% --------
% canExchange

deprecationWarning('canProduce','canExchange(model,''produce'',...)');

produced=canExchange(model,'produce',varargin{:});
end
