function followChanged(model,fluxesA,fluxesB, cutOffChange, cutOffFlux, cutOffDiff, metaboliteList)
% followChanged
%	Prints fluxes and reactions for each of the reactions that results in
%   different fluxes compared to the reference case.
%
%   model           a model structure
%   fluxesA         flux vector for the test case
%   fluxesB         flux vector for the reference test
%   cutOffChange	reactions where the fluxes differ by less than
%                   this many percent won't be printed (opt, default 10^-8)
%   cutOffFlux      reactions where the absolute value of both fluxes
%                   are below this value won't be printed (opt,
%                   default 10^-8)
%   cutOffDiff      reactions where the fluxes differ by less than
%                   cutOffDiff won't be printed (opt, default 10^-8)
%   metaboliteList  cell array of metabolite names. Only reactions
%                   involving any of these metabolites will be
%                   printed (opt)
%
%   Usage: followChanged(model,fluxesA,fluxesB, cutOffChange, cutOffFlux,
%           cutOffDiff, metaboliteList)

%Checks if a cut off flux has been set
if nargin<4
    cutOffChange=10^-8;
end
if nargin<5
    cutOffFlux=10^-8;
end
if nargin<6
    cutOffDiff=10^-8;
end
if nargin<7
    metaboliteList=[];
else
    metaboliteList=convertCharArray(metaboliteList);
end

%If a metabolite list is to be used, then find all the reactions involving
%any of those metabolites Finds the metabolites
if nargin>6
    reactionIndexes=[];
    for i=1:length(metaboliteList)
        metaboliteIndex=find(strcmpi(metaboliteList(i),model.metNames)); %Should use id maybe, setting
        if ~isempty(metaboliteIndex)
            [~, b]=find(model.S(metaboliteIndex,:));
            reactionIndexes=[reactionIndexes; b(:)];
        else
            fprintf('Could not find any reactions with the metabolite %s\n\n',char(metaboliteList(i)))
        end
    end
    reactionIndexes=unique(reactionIndexes);
else
    reactionIndexes=(1:length(fluxesA))';
end

%Finds the reactions where either flux is at or above the cutOffFlux value
in1=find(abs(fluxesA(reactionIndexes))>=cutOffFlux)';
in2=find(abs(fluxesB(reactionIndexes))>=cutOffFlux)';
ineither=reactionIndexes(unique([in1 in2]));

%Keep only those solutions where the difference is larger than or equal to
%cutOffDiff
ineither=ineither(find(abs(fluxesA(ineither)-fluxesB(ineither))>=cutOffDiff));

%Finds the reactions where the fluxes differ more than cutOffChange percent
%First check those fluxes that are non-zero in solution1.x
nonZeroFluxes=ineither(find(fluxesA(ineither)));
quota=1+cutOffChange/100;
larger=nonZeroFluxes(find((fluxesB(nonZeroFluxes)./fluxesA(nonZeroFluxes))>=(quota)))';
smaller=nonZeroFluxes(find((fluxesB(nonZeroFluxes)./fluxesA(nonZeroFluxes))<(1/quota)))';
fluxIndexes=[larger smaller];

%Then add those where solution1 has a zero flux
zeroFluxes=ineither(find(fluxesA(ineither)==0));
fluxIndexes=unique([fluxIndexes zeroFluxes(find(abs(fluxesB(zeroFluxes))>=cutOffFlux))']);

formulas=constructEquations(model,model.rxns(fluxIndexes));

if nargin>4
    if nargin>5
        fprintf('These reactions have flux values that differ by more than %s percent, absolute values above %s, and a total difference above %s (%s reactions)\n\n',num2str(cutOffChange),num2str(cutOffFlux),num2str(cutOffDiff),num2str(length(formulas)));
    else
        fprintf('These reactions have flux values that differ by more than %s percent and absolute values above %s (%s reactions)\n\n',num2str(cutOffChange),num2str(cutOffFlux),num2str(length(formulas)));
    end
else
    fprintf('These reactions have flux values that differ by more than %s percent (%s reactions)\n\n',num2str(cutOffChange),num2str(length(formulas)));
end

metaboliteNames=[];
for i=1:length(metaboliteList)
    metaboliteNames=[metaboliteNames char(metaboliteList(i)) ' '];
end

if ~isempty(metaboliteNames)
    fprintf('Only prints reactions involving one or more of the following metabolites:\n%s\n\n',metaboliteNames)
end

for i=1:length(formulas)
    fluxText=['Flux: ' num2str(fluxesA(fluxIndexes(i))) ' Reference flux: ' num2str(fluxesB(fluxIndexes(i))) ' Difference: ' num2str(fluxesA(fluxIndexes(i))-fluxesB(fluxIndexes(i)))];
    fprintf('%s: %s\n\t%s\n\t%s\n\n', char(model.rxns(fluxIndexes(i))), char(formulas(i)),...
        char(model.rxnNames(fluxIndexes(i))),fluxText);
end
end
