function gprob = cobraToGurobiProb(prob)

gprob.obj = prob.c;
gprob.rhs = prob.b;
gprob.lb = prob.lb;
gprob.ub = prob.ub;
gprob.A = prob.A;
gprob.sense = '=';
gprob.modelsense = 'min';
gprob.vtype = repmat('C', 1, size(prob.A, 2));
if isfield(prob,'basis') && ~isempty(prob.basis)
    gprob.cbasis=prob.basis.cbasis;
    gprob.vbasis=prob.basis.vbasis;
end
end