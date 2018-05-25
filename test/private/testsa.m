function ST=testsa(S)
% N=1e4;
% S=sparse(rand(N,N));
% T=sparse(rand(N,N));

% ts2=tic;
% asd2=-sparse(eye(size(S,1)));
% te2=toc(ts2)
% ts1=tic;
% asd=sparse([S T]);
% te1=toc(ts1)


% ts3=tic;
% [is, js, ss] = find(S);
% [it, jt, st] = find(T);
% ST = sparse([is; it + size(S,1)], [js; jt], [ss; st]);
% te3=toc(ts3)
T=-eye(size(S,1));

tic
for i=1:1
	[is, js, ss] = find(S);
	[it, jt, st] = find(T);
	ST = sparse([is;it], [js;jt+size(S,2)], [ss;st]);
end
toc

tic

for i=1:1 [S T]; end
toc
end
