%These tests are using the Human-GEM model.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% TC0901 - Test that all reactions in the minimal model to be
%          used by the MILP can carry flux > 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%

%In the MILP, the flux of reactions with negative rxn scores cannot exceed 100
%At the same time, the flux of positive reactions need to be >= 0.1 to be able
%to turn them on. Likewise, we force the flux on essential reactions to be >= 0.1.
%Therefore, it is important that all reactions can get a flux > 0.1, while all
%reactions are limited to a flux <= 100. We can test this in a similar way to the
%haveFlux function in RAVEN, so that code is copied here and modified.

%This code is partly copied from the function haveFlux
%The prepData takes a couple of hours to generate. If you have it somewhere,
%better to load it
%cd C:\Work\MatlabCode\components\human-GEM\Human-GEMftINIT\Human-GEM %replace this with the root of the Human-GEM repo.
%ihuman = importYaml('model/Human-GEM.yml');
%save('model/Human-GEM.mat', 'ihuman');
%load('model/Human-GEM.mat')
%prepData = prepHumanModelForftINIT(ihuman, false);
prepData2 = load('model/PrepData2.mat').prepData;

upperLimit = 100;
lowerLimit = 0.1;

model = prepData2.minModel;
model = prepData.minModel;




model.lb(model.lb==-inf)=-upperLimit;
model.ub(model.ub==inf)=upperLimit;
model.lb(model.lb==-1000)=-upperLimit;
model.ub(model.ub==1000)=upperLimit;

%check
unique(model.lb) %ok, only 0 and -100
unique(model.ub) %ok, only 100

%First make a loop where we optimize for all
onFwd = false(numel(model.rxns,1));
onRev = false(numel(model.rxns,1));

iter = 1;
while true
    disp([num2str(iter) ': ' num2str(sum(onFwd))])
    model.c=ones(numel(model.c),1);
    model.c(onFwd) = 0; %don't include the ones already verified
    
    sol=solveLP(model);
    if isempty(sol.x)
        disp('Failed');
        break;
    end
    
    currRes = sol.x > lowerLimit;
    if sum(currRes & ~onFwd) == 0
        break;%we didn't turn any more rxns on
    end
    
    onFwd = onFwd | currRes;
    iter = iter + 1;
end

sum(~onFwd)%zero, perfect

%Now the other direction (reversible):
onRev = false(numel(model.rxns,1));

iter = 1;
while true
    disp([num2str(iter) ': ' num2str(sum(onRev))])
    model.c=-ones(numel(model.c),1);
    model.c(~model.rev) = 0;
    model.c(onRev) = 0; %don't include the ones already verified
    
    sol=solveLP(model);
    if isempty(sol.x)
        disp('Failed');
        break;
    end
    
    currRes = sol.x < -lowerLimit;
    if sum(currRes & ~onRev) == 0
        break;%we didn't turn any more rxns on
    end
    
    onRev = onRev | currRes;
    iter = iter + 1;
end

sum(~onRev(model.rev == 1))%

onRev = false(numel(model.rxns,1));

indLeft = find(model.rev == 1 & ~onRev);

constructEquations(model, model.rxns(indLeft))%no scaling here...

%go trough the last one by one
for i = 1:numel(indLeft)
    model.c=zeros(numel(model.c),1);
    model.c(indLeft(i)) = -1;
    
    sol=solveLP(model);
    if isempty(sol.x)
        disp('Failed');
        break;
    end
    
    currRes = sol.x(indLeft(i)) < -lowerLimit;
    disp([num2str(indLeft(i)) ': ' num2str(currRes)])
    onRev(indLeft(i)) = currRes;
    onRev = onRev | currRes;
    iter = iter + 1;
end

sum(~onRev(model.rev == 1))%0, perfect. This means that we can use the limits 0.1 and 100 in the MILP!


