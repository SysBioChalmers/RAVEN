function outProb = mosekToCobraProb(inProb)
	outProb = [];

	outProb.c = [inProb.c;zeros(size(inProb.a,1),1)];
	outProb.A = [inProb.a -speye(size(inProb.a,1))]; 
	outProb.b = zeros(size(inProb.a,1), 1);
	outProb.lb = [inProb.blx; inProb.blc]; 
	outProb.ub = [inProb.bux; inProb.buc];
	outProb.osense = 1;
	outProb.csense = char(zeros(size(inProb.a,1),1));
	outProb.csense(:) = 'E';

	% milp 
	if(isfield(inProb,'ints'))
		outProb.vartype = repmat('C', 1, size(outProb.A, 2)); 
		outProb.vartype(inProb.ints.sub) = 'B';
		outProb.x0=[]; 
	end

end
