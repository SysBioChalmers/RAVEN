function errorFlag=followFluxes(model, fluxesA, lowerFlux, upperFlux, fluxesB)
% followFluxes
%	Prints fluxes and reactions for each of the reactions that results in
%   fluxes in the specified interval.
%
%   model       a model structure
%   fluxesA     flux vector for the test case
%   lowerFlux	only reactions with fluxes above this cutoff
%               value are displayed
%   upperFlux   only reactions with fluxes below this cutoff
%               value are displayed (opt, default Inf)
%   fluxesB     flux vector for the reference case(opt)
%
%   Usage: errorFlag=followFluxes(model, fluxesA, lowerFlux, upperFlux,
%           fluxesB)
%
%   Rasmus Agren, 2010-12-16
%

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