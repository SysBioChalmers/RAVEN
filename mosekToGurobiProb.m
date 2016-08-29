function gprob = mosekToGurobiProb(prob)
	gprob.obj=[prob.c;zeros(size(prob.a,1),1)];
	gprob.A = [prob.a -eye(size(prob.a,1))];
	gprob.rhs = zeros(size(prob.a,1), 1);
	gprob.lb = [prob.blx; prob.blc];
	gprob.ub = [prob.bux; prob.buc];
	gprob.sense = '=';
	gprob.modelsense = 'min';
	gprob.vtype = repmat('C', 1, size(gprob.A, 2)); 

	% the binary type variables must be defined for milp 
	if(isfield(prob,'ints')) gprob.vtype(prob.ints.sub) = 'B'; end

	% hotstart
	if(isfield(prob,'sol'))
		if(isfield(prob.sol,'bas'))
			gprob.vbasis=prob.sol.bas.vbasis;
			gprob.cbasis=prob.sol.bas.cbasis;
		end
	end
end