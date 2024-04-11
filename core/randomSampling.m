function [solutions, goodRxns]=randomSampling(model,nSamples,replaceBoundsWithInf,supressErrors,runParallel,goodRxns,minFlux)
% randomSampling
%   Performs random sampling of the solution space, as described in Bordel
%   et al. (2010) PLOS Compt Biol (doi:10.1371/journal.pcbi.1000859).
%
% Input:
%   model                   a model structure
%   nSamples                the number of solutions to return
%                           (opt, default 1000)
%   replaceBoundsWithInf    replace the largest upper bounds with Inf and
%                           the smallest lower bounds with -Inf. This is
%                           needed in order to get solutions without loops
%                           if your model has for example 1000/-1000 as
%                           arbitary large bounds. If your model only has
%                           "biologically relevant" bounds, then set this
%                           to false (opt, default true)
%   supressErrors           the program will halt if it has problems
%                           finding non-zero solutions which are not
%                           involved in loops. This could be because the
%                           constraints on the model are too relaxed (such
%                           as unlimited glucose uptake) or too strict
%                           (such as too many and too narrow constraints)
%                           (opt, default false)
%   runParallel             speed up calculations by parallel processing.
%                           Requires MATLAB Parallel Computing Toolbox. If
%                           this is not installed, the calculations will
%                           not be parallelized, regardless what is
%                           indicated as runParallel. (opt, default true)
%   goodRxns                double vector of indexes of those reactions
%                           that are not involved in loops and can be used
%                           as random objective functions, as generated by
%                           a previous run of randomSampling on the same
%                           model (opt, default empty)
%   minFlux                 determines if a second optimization should be
%                           performed for each random sample, to minimize
%                           the number of fluxes and thereby preventing
%                           loops. Typically, loops are averaged out when a
%                           large number of samples are taken, but this is
%                           not always the case (opt, default false)
%
%
% Output:
%   solutions               matrix with the solutions
%   goodRxns                double vector of indexes of those reactions
%                           that are not involved in loops and can be used
%                           as random objective functions
%
%   The solutions are generated by maximizing (with random weights) for a
%   random set of three reactions. For reversible reactions it randomly
%   chooses between maximizing and minimizing.
%
%   Usage: solutions=randomSampling(model,nSamples,replaceBoundsWithInf,supressErrors,runParallel,goodRxns,minFlux)

if nargin<2 | isempty(nSamples)
    nSamples=1000;
end
if nargin<3 | isempty(replaceBoundsWithInf)
    replaceBoundsWithInf=true;
end
if nargin<4 | isempty(supressErrors)
    supressErrors=false;
end
if nargin<5 | isempty(runParallel)
    runParallel=true;
end
if nargin<7 | isempty(minFlux)
    minFlux=false;
end

[ps, oldPoolAutoCreateSetting] = parallelPoolRAVEN(runParallel);

nObjRxns=2; %Number of reactions in the objective function in each iteration

%Simplify the model to speed stuff up a little. Keep original mapping
originalRxns=model.rxns;
model=simplifyModel(model,false,false,true,true);

%Then change the bounds to +/- Inf. This is needed in order to not have
%loops in the solutions
if replaceBoundsWithInf==true
    model.ub(model.ub==max(model.ub))=Inf;
    if min(model.lb)<0 % Only negative lower bounds should be set to -Inf
        model.lb(model.lb==min(model.lb))=-Inf;
    end
end

%Check that the model is feasible given the constraints
[sol,~]=solveLP(model);
if isempty(sol.x)
    EM='The model has no feasible solution, likely due to incompatible constraints';
    dispEM(EM);
elseif sol.f==0 && showProgress
    warning('The model objective function cannot reach a non-zero value. This might be intended, so randomSampling will continue, but this could indicate problems with your model')
end

[~,hsSol]=solveLP(model);
nRxns = numel(model.rxns);
%Reactions which can be involved in loops should not be optimized for.
%Check which reactions reach an arbitary high upper bound
if nargin<6 || isempty(goodRxns)
    goodRxns = true(numel(model.rxns),numel(model.rxns));
    goodRxns = num2cell(goodRxns,1);
    PB = ProgressBar2(nRxns,'Prepare goodRxns not involved in loops','cli');
    parfor (i=1:nRxns)
        testModel=setParam(model,'eq',model.rxns(i),1000);
        sol=solveLP(testModel,0,[],hsSol);
        if ~isempty(sol.f)
            notGood = abs(sol.x)>999;
            goodRxns{i}(notGood)=false;
        else
            %If the reaction is reversible, also check in that direction
            if model.rev(i)
                testModel=setParam(model,'eq',model.rxns(i),-1000);
                sol=solveLP(testModel,0,[],hsSol);
                if ~isempty(sol.f)
                    notGood = abs(sol.x)>999;
                    goodRxns{i}(notGood)=false;
                end
            end
        end
        count(PB);
    end
    goodRxns = cell2mat(goodRxns);
    goodRxns = find(prod(goodRxns,2));
end

%Reserve space for a solution matrix
sols = zeros(numel(model.rxns),nSamples);
sols = num2cell(sols,1);

%Main loop
PB = ProgressBar2(nSamples,'Performing random sampling','cli');
parfor i=1:nSamples
    badSolutions = 1;
    tmpModel = model;
    while lt(0,badSolutions)
        rxns = randsample(numel(goodRxns),nObjRxns);
        tmpModel.c = zeros(numel(tmpModel.rxns),1);
        multipliers = randsample([-1 1],nObjRxns,true);
        multipliers(tmpModel.rev(goodRxns(rxns))==0) = 1;
        tmpModel.c(goodRxns(rxns)) = rand(nObjRxns,1).*multipliers;
        if true(minFlux)
            sol=solveLP(tmpModel,1,[],hsSol);
        else
            sol=solveLP(tmpModel,0,[],hsSol);
        end
        if any(sol.x) && abs(sol.f)>10^-8
            sols{i} = sol.x;
            badSolutions = 0;
        else
            badSolutions = badSolutions+1;
            %If it only finds bad solutions then throw an error.
            if badSolutions == 100 && supressErrors==false
                error('The program is having problems finding non-zero solutions (ignoring reactions that might be involved in loops). Review the constraints on your model. Set supressErrors to true to ignore this error');
            end
        end
    end
    count(PB);
end

%Map to original model
sols = cell2mat(sols);
[~, I]=ismember(model.rxns,originalRxns);
solutions=zeros(numel(originalRxns),nSamples);
solutions(I,:)=sols;
solutions=sparse(solutions);
% Reset original Parallel setting
ps.Pool.AutoCreate = oldPoolAutoCreateSetting;
end

%To use instead of the normal Matlab randsample function. This is in order
%to not depend on the Matlab statistical toolbox.
function I=randsample(n,k,replacement)
if nargin<3
    replacement=false;
end
%n can be a integer, which leads to I being sampled from 1:n, or it can be
%a population to sample from.
if numel(n)==1 && isnumeric(n)
    n=1:n;
end
%Loop and get random numbers until the list is unique. This is only a good
%option is the number of samples is small compared to the population. There
%are several checks that should be made here, for example regarding size
%and that the number of samples is <=population size if replacement==false.
%This is not the case in randomSampling, so such checks are ignored
while true
    J=randi(numel(n),[k,1]);
    if replacement==true || numel(J)==numel(unique(J))
        I=n(J);
        break;
    end
end
I=I(:);
end
