function errorFlag=followFluxes(model, fluxesA, lowerFlux, upperFlux, fluxesB)
% followFluxes  Print reactions with fluxes in a specified interval.
%
% Prints fluxes and reactions for each of the reactions that result in
% fluxes within the specified interval.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% fluxesA : double
%     flux vector for the test case.
% lowerFlux : double
%     only reactions with fluxes above this cutoff value are displayed.
% upperFlux : double, optional
%     only reactions with fluxes below this cutoff value are displayed
%     (default Inf).
% fluxesB : double, optional
%     flux vector for the reference case.
%
% Returns
% -------
% errorFlag : double
%     set to 1 if upperFlux is not larger than lowerFlux, otherwise empty.
%
% Examples
% --------
%     errorFlag=followFluxes(model,fluxesA,lowerFlux,upperFlux,fluxesB);

%Checks that the upper flux is larger than the lower flux
if nargin>3
    if upperFlux<=lowerFlux
        errorFlag=1;
        return;
    end
end

%Gets the fluxes for the reactions
if nargin<4
    fluxIndexes=find(fluxesA>=lowerFlux);
else
    fluxIndexes=find(fluxesA>=lowerFlux & fluxesA<=upperFlux);
end

%Finds the involved reactions
formulas = constructEquations(model,model.rxns(fluxIndexes));

if nargin>3
    fprintf('These reactions have flux values between %s and %s\n\n',num2str(lowerFlux),num2str(upperFlux));
else
    fprintf('These reactions have flux values above %s\n\n',num2str(lowerFlux));
end
for i=1:length(formulas)
    if nargin>4
        fluxText=['Flux: ' num2str(fluxesA(fluxIndexes(i))) ' Reference flux: ' num2str(fluxesB(fluxIndexes(i)))];
    else
        fluxText=['Flux: ' num2str(fluxesA(fluxIndexes(i)))];
    end
    fprintf('%s: %s\n\t%s\n\t%s\n', char(model.rxns(fluxIndexes(i))), char(formulas(i)),...
        char(model.rxnNames(fluxIndexes(i))),fluxText);
end
end
