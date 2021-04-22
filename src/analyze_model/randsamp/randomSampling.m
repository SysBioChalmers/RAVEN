function [solutions, goodRxns]=randomSampling(model,nSamples,replaceBoundsWithInf,supressErrors,showProgress,goodRxns,minFlux)
% randomSampling
%   Returns a number of random solutions
%
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
%   showProgress            if true, it will display in the command window 
%                           how many iterations have been done (opt, default
%                           false)
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
%   solutions               matrix with the solutions
%   goodRxns                double vector of indexes of those reactions
%                           that are not involved in loops and can be used
%                           as random objective functions
%
%   The solutions are generated by maximizing (with random weights) for a
%   random set of three reactions. For reversible reactions it randomly
%   chooses between maximizing and minimizing.
%
%   Usage: solutions=randomSampling(model,nSamples,replaceBoundsWithInf,supressErrors,showProgress,goodRxns,minFlux)

if nargin<2
    nSamples=1000;
end
if nargin<3
    replaceBoundsWithInf=true;
end
if nargin<4
    supressErrors=false;
end
if nargin<5
    showProgress=false;
end
if nargin<6
    minFlux=false;
end

nRxns=2; %Number of reactions in the objective function in each iteration

%First check that the model is feasible given the constraints
sol=solveLP(model);
if isempty(sol.x)
    EM='The model has no feasible solution';
    dispEM(EM);
end

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

%Reactions which can be involved in loops should not be optimized for.
%Check which reactions reach an arbitary high upper bound
if nargin<6 || isempty(goodRxns)
    goodRxns=true(numel(model.rxns),1);
    for i=1:numel(model.rxns)
        if showProgress && rem(i,100) == 0
            disp(['Preparing random sampling: ready with ' num2str(i) '/' num2str(numel(model.rxns)) ' rxns'])
        end
        if goodRxns(i)==true
            testModel=setParam(model,'eq',model.rxns(i),1000);
            sol=solveLP(testModel);
            if ~isempty(sol.f)
                goodRxns(abs(sol.x)>999)=false;
            else
                %If the reaction is reversible, also check in that direction
                if model.rev(i)
                    testModel=setParam(model,'eq',model.rxns(i),-1000);
                    sol=solveLP(testModel);
                    if ~isempty(sol.f)
                        goodRxns(abs(sol.x)>999)=false;
                    end
                end
            end
        end
    end
    goodRxns=find(goodRxns);
end

%Reserve space for a solution matrix
sols=zeros(numel(model.rxns),nSamples);

%Main loop
counter=1;
badSolutions=0;
while counter<=nSamples
    if showProgress && rem(counter,100) == 0
        disp(['Performing random sampling: ready with ' num2str(counter) '/' num2str(nSamples) ' iterations'])
    end
    rxns=randsample(numel(goodRxns),nRxns);
    model.c=zeros(numel(model.rxns),1);
    multipliers=randsample([-1 1],nRxns,true);
    multipliers(model.rev(goodRxns(rxns))==0)=1;
    model.c(goodRxns(rxns))=rand(nRxns,1).*multipliers;
    if true(minFlux)
        sol=solveLP(model,1);
    else
        sol=solveLP(model);
    end
    if any(sol.x)
        if abs(sol.f)>10^-8
            sols(:,counter)=sol.x;
            counter=counter+1;
            badSolutions=0;
        else
            badSolutions=badSolutions+1;
            %If it only finds bad solutions then throw an error.
            if badSolutions==50 && supressErrors==false
                EM='The program is having problems finding non-zero solutions that are not involved in loops. Review the constraints on your model. Set supressErrors to true to ignore this error';
                dispEM(EM);
            end
        end
    end
end

%Map to original model
[~, I]=ismember(model.rxns,originalRxns);
solutions=zeros(numel(originalRxns),nSamples);
solutions(I,:)=sols;
solutions=sparse(solutions);
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
