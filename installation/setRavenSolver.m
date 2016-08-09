function out = setRavenSolver(solver)
	if (~ischar(solver)) dispEM('Input should be a string.'); end
	global RAVENSOLVER;
	RAVENSOLVER=solver;
end