function gprob = cobraToGurobiProb(prob)

gprob.obj = prob.c;
gprob.rhs = prob.b;
gprob.lb = prob.lb;
gprob.ub = prob.ub;
gprob.A = prob.A;
gprob.sense = '=';
gprob.modelsense = 'min';

% the binary type variables must be defined for milp
% use this more complicated approach to ignore "vartype" capitalization
f = fieldnames(prob);
[~,vartype_index] = ismember('vartype', lower(f));
if vartype_index > 0
    gprob.vtype = prob.(f{vartype_index});
else
    gprob.vtype = repmat('C', 1, size(prob.A, 2));
end

if isfield(prob,'basis') && ~isempty(prob.basis)
    gprob.cbasis=prob.basis.cbasis;
    gprob.vbasis=prob.basis.vbasis;
end
end