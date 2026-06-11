function report=checkRxn(model,rxn,varargin)
% checkRxn  Check which reactants can be synthesized and products consumed.
%
% Checks which reactants in a reaction that can be synthesized and which
% products that can be consumed. This is primarily for debugging reactions
% which cannot have flux.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxn : char
%     the id of one reaction to check.
%
% Name-Value Arguments
% --------------------
% cutoff : double
%     minimal flux for successful production/consumption (default 10^-7).
% revDir : logical
%     true if the reaction should be reversed (default false).
% printReport : logical
%     print a report (default true).
%
% Returns
% -------
% report : struct
%     report with fields:
%
%     - reactants : array with reactant indexes
%     - canMake : boolean array, true if the corresponding reactant can be
%       synthesized by the rest of the metabolic network
%     - products : array with product indexes
%     - canConsume : boolean array, true if the corresponding product can be
%       consumed by the rest of the metabolic network
%
% Examples
% --------
%     report = checkRxn(model, rxn, cutoff, revDir, printReport);

rxn=char(rxn);
p=parseRAVENargs(varargin, {'cutoff',[]; ...
    'revDir',false; ...
    'printReport',true});
cutoff=p.cutoff;
revDir=p.revDir;
printReport=p.printReport;

if isempty(cutoff)
    cutoff=10^-7;
end

[I, rxnID]=ismember(rxn,model.rxns);

if ~I
    EM='Reaction ID not found';
    dispEM(EM);
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
    [tempModel, testRxn]=addExchangeRxns(model,'out',report.reactants(i));
    tempModel=setParam(tempModel,'obj',testRxn,1);
    sol=solveLP(tempModel);
    if sol.f>cutoff
        report.canMake(i)=true;
    else
        if printReport==true
            fprintf(['Failed to make ' model.metNames{report.reactants(i)} '[' model.comps{model.metComps(report.reactants(i))} ']\n']);
        end
    end
end

for i=1:numel(report.products)
    [tempModel, testRxn]=addExchangeRxns(model,'in',report.products(i));
    tempModel=setParam(tempModel,'obj',testRxn,1);
    sol=solveLP(tempModel);
    if sol.f>cutoff
        report.canConsume(i)=true;
    else
        if printReport==true
            fprintf(['Failed to consume ' model.metNames{report.products(i)} '[' model.comps{model.metComps(report.products(i))} ']\n']);
        end
    end
end
end
