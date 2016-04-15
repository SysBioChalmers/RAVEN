function report=checkRxn(model,rxn,cutoff,revDir,printReport)
% checkRxn
%   Checks which reactants in a reaction that can be synthesized and which
%   products that can be consumed. This is primarily for debugging
%   reactions which cannot have flux
%
%   model       a model structure
%   rxn         the id of the reaction to check
%   cutoff      minimal flux for successful production/consumption (opt,
%               default 10^-7) 
%   revDir      true if the reaction should be reversed (opt, default
%               false)
%   printReport print a report (opt, default true)
%
%   report
%       reactants   array with reactant indexes
%       canMake     boolean array, true if the corresponding reactant can be
%                   synthesized
%       products    array with product indexes
%       canConsume  boolean array, true if the corresponding reactant can
%                   be consumed
%
%   Usage: report=checkRxn(model,rxn,cutoff,revDir,printReport)
%
%   Rasmus Agren, 2013-09-17
%

%Convert to cell string
if isstr(rxn)
    rxn={rxn};
end
if nargin<3
    cutoff=10^-7;
end
if nargin<4
    revDir=false;
end
if isempty(cutoff)
    cutoff=10^-7;
end
if nargin<5
    printReport=true;
end

[I rxnID]=ismember(rxn,model.rxns);

if ~I
    dispEM('Reaction ID not found');
end

if revDir==false
    report.reactants=find(model.S(:,rxnID)<0);
    report.products=find(model.S(:,rxnID)>0);
else
    report.reactants=find(model.S(:,rxnID)>0);
    report.products=find(model.S(:,rxnID)<0);
end
report.canMake=false(numel(report.reactants),1);
report.canConsume=false(numel(report.products),1);

%Remove this field as it would give an annoying note otherwise
if isfield(model,'rxnComps')
    model=rmfield(model,'rxnComps');
end

%There are several ways to do this. Here I choose to add the reactions one
%by one and checking their bounds. This might not be optimal
for i=1:numel(report.reactants)
    [tempModel testRxn]=addExchangeRxns(model,'out',report.reactants(i));
    tempModel=setParam(tempModel,'obj',testRxn,1);
    sol=solveLP(tempModel);
    if sol.f*-1>cutoff
       report.canMake(i)=true; 
    else
        if printReport==true
           fprintf(['Failed to make ' model.metNames{report.reactants(i)} '[' model.comps{model.metComps(report.reactants(i))} ']\n']);
        end
    end
end

for i=1:numel(report.products)
    [tempModel testRxn]=addExchangeRxns(model,'in',report.products(i));
    tempModel=setParam(tempModel,'obj',testRxn,1);
    sol=solveLP(tempModel);
    if sol.f*-1>cutoff
       report.canConsume(i)=true; 
    else
        if printReport==true
            fprintf(['Failed to consume ' model.metNames{report.products(i)} '[' model.comps{model.metComps(report.products(i))} ']\n']);
        end
    end
end
end
