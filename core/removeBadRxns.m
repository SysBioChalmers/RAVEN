function [newModel, removedRxns]=removeBadRxns(model,rxnRules,ignoreMets,isNames,balanceElements,refModel,ignoreIntBounds,printReport)
% removeBadRxns
%   Iteratively removes reactions which enable production/consumption of some
%   metabolite without any uptake/excretion
%
%   model                   a model structure. For the intented function,
%                           the model shouldn't allow for any uptake/excretion.
%                           The easiest way to achieve this is to import the
%                           model using importModel('filename',false)
%   rxnRules                1: only remove reactions which are unbalanced
%                           2: also remove reactions which couldn't be checked for
%                           mass balancing
%                           3: all reactions can be removed
%                           (opt, default 1)
%   ignoreMets              either a cell array of metabolite IDs, a logical vector
%                           with the same number of elements as metabolites in the model,
%                           of a vector of indexes for metabolites to exclude from
%                           this analysis (opt, default [])
%   isNames                 true if the supplied mets represent metabolite names
%                           (as opposed to IDs). This is a way to delete
%                           metabolites in several compartments at once without
%                           knowing the exact IDs. This only works if ignoreMets
%                           is a cell array (opt, default false)
%   balanceElements         a cell array with the elements for which to
%                           balance the reactions. May contain any
%                           combination of the elements defined in parseFormulas
%                           (opt, default {'C';'P';'S';'N';'O'})
%   refModel                a reference model which can be used to ensure
%                           that the resulting model is still functional.
%                           The intended use is that the reference model is
%                           a copy of model, but with uptake/excretion allowed and
%                           some objectives (such as production of biomass)
%                           constrained to a non-zero flux. Before a
%                           reaction is removed from "model" the function first
%                           checks that the same deletion in "refModel"
%                           doesn't render the problem unfeasible (opt)
%   ignoreIntBounds         true if internal bounds (including reversibility)
%                           should be ignored. Exchange reactions are not affected.
%                           This can be used to find unbalanced solutions which are
%                           not possible using the default constraints (opt,
%                           default false)
%   printReport             true if a report should be printed (opt,
%                           default false)
%
%   newModel               	a model structure after the problematic
%                           reactions have been deleted
%   removedRxns             a cell array with the reactions that were
%                           removed
%
%   The purpose of this function is to remove reactions which enable
%   production/consumption of metabolites even when exchange reactions aren't used.
%   Many models, especially if they are automatically inferred from
%   databases, will have unbalanced reactions which allow for
%   net-production/consumption of metabolites without any consumption/excretion.
%   A common reason for this is when general compounds have different meaning
%   in different reactions (as DNA has in these two reactions).
%       dATP + dGTP + dCTP + dTTP <=> DNA + 4 PPi
%       0.25 dATP + 0.25 dGTP + 0.25 dCTP + 0.25 dTTP <=> DNA + PPi
%   Reactions that are problematic like this are always elementally
%   unbalanced, but it is not always that you would like to exclude all
%   unbalanced reactions from your model.
%   This function tries to remove as few problematic reactions as possible
%   so that the model cannot produce/consume anything from nothing. This is done by
%   repeatedly calling makeSomething/consumeSomething, checking if any of
%   the involved reactions are elementally unbalanced, remove one of them,
%   and then iterating until no metabolites can be produced/consumed.
%   makeSomething is called before consumeSomething.
%
%   Usage: [newModel, removedRxns]=removeBadRxns(model,rxnRules,...
%       ignoreMets,isNames,refModel,ignoreIntBounds,printReport)

if nargin<2
    rxnRules=1;
end
if nargin<3
    ignoreMets=[];
elseif ~islogical(ignoreMets) && ~isnumeric(ignoreMets)
    ignoreMets=convertCharArray(ignoreMets);
end
if nargin<4
    isNames=false;
end
if nargin<5
    balanceElements={'C';'P';'S';'N';'O'};
else
    balanceElements=convertCharArray(balanceElements);
end
if nargin<6
    refModel=[];
else
    if ~isempty(refModel)
        EM='The feature to supply a reference model is currently not supported';
        dispEM(EM,false);
    end
end
if nargin<7
    ignoreIntBounds=false;
end
if nargin<8
    printReport=false;
end

%Check if the model has open exchange reactions and print a warning in that
%case
if ~isfield(model,'unconstrained')
    [~, I]=getExchangeRxns(model);
    if any(I)
        if any(model.lb(I)~=0) || any(model.ub(I)~=0)
            EM='The model contains open exchange reactions. This is not the intended use of this function. Consider importing your model using importModel(filename,false)';
            dispEM(EM,false);
        end
    end
end

%Check that the model is feasible
sol=solveLP(model);
if isempty(sol.f)
    EM='The model is not feasible. Consider removing lower bounds (such as ATP maintenance)';
    dispEM(EM);
end

%Check that the reference model is feasible
if any(refModel)
    sol=solveLP(refModel);
    if isempty(sol.f)
        EM='The reference model is not feasible';
        dispEM(EM);
    end
end

%Initialize some stuff
removedRxns={};
balanceStructure=getElementalBalance(model);

%Check which elements to balance for
if ~isempty(setdiff(balanceElements,balanceStructure.elements.abbrevs))
    EM='Could not recognize all elements to balance for';
    dispEM(EM);
end
bal=ismember(balanceStructure.elements.abbrevs,balanceElements);

%Main loop. First run for makeSomething, second for consumeSomething
warned=false(2,1); %This is to prevent the same warning being printed multiple times if rxnRules==3
for i=1:2
    while 1
        %Make some metabolite using as few reactions as possible
        if i==1
            [solution, metabolite]=makeSomething(model,ignoreMets,isNames,false,true,[],ignoreIntBounds);
            if ~isempty(solution)
                if printReport
                    fprintf(['Can make: ' model.metNames{metabolite(1)} '\n']);
                end
            else
                %If no solution could be found, then finish
                break;
            end
        else
            [solution, metabolite]=consumeSomething(model,ignoreMets,isNames,false,[],ignoreIntBounds);
            if ~isempty(solution)
                if printReport
                    fprintf(['Can consume: ' model.metNames{metabolite(1)} '\n']);
                end
            else
                %If no solution could be found, then finish
                break;
            end
        end
        
        %Find all reactions that are unbalanced and still carry flux
        I=find(abs(solution)>10^-8 & balanceStructure.balanceStatus>=0 & ~all(balanceStructure.leftComp(:,bal)==balanceStructure.rightComp(:,bal),2));
        
        %If there are unbalanced rxns then delete one of them and iterate
        if any(I)
            rxnToRemove=I(randsample(numel(I),1));
        else
            %If there are no unbalanced rxns in the solution
            if rxnRules==1
                if i==1
                    EM=['No unbalanced reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 2 or 3 for a more exhaustive search'];
                    dispEM(EM,false);
                else
                    EM=['No unbalanced reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 2 or 3 for a more exhaustive search'];
                    dispEM(EM,false);
                end
                break;
            else
                %Find reactions which are not checked for mass balancing,
                %but that still carry flux
                I=find(abs(solution)>10^-8 & balanceStructure.balanceStatus<0);
                
                %If there are any such reactions, remove one of them and
                %iterate
                if any(I)
                    rxnToRemove=I(randsample(numel(I),1));
                else
                    if rxnRules==2
                        %This happens when all reactions used are balanced
                        %according to the metabolite formulas. This cannot
                        %be, and indicates that one or more of the formulas
                        %are wrong. Print a warning and delete any reaction
                        %with flux
                        if i==1
                            EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 3 for a more exhaustive search'];
                            dispEM(EM,false);
                        else
                            EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 3 for a more exhaustive search'];
                            dispEM(EM,false);
                        end
                        break;
                    else
                        if i==1
                            if warned(1)==false
                                EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". This indicates some error in the metabolite formulas. Removing random reactions in the solution'];
                                dispEM(EM,false);
                                warned(1)=true;
                            end
                        else
                            if warned(2)==false
                                EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". This indicates some error in the metabolite formulas. Removing random reactions in the solution'];
                                dispEM(EM,false);
                                warned(2)=true;
                            end
                        end
                        I=find(abs(solution)>10^-8);
                        rxnToRemove=I(randsample(numel(I),1));
                    end
                end
            end
        end
        removedRxns=[removedRxns;model.rxns(rxnToRemove)];
        if printReport
            fprintf(['\tRemoved: '  model.rxns{rxnToRemove} '\n']);
        end
        model=removeReactions(model,rxnToRemove);
        balanceStructure.balanceStatus(rxnToRemove)=[];
        balanceStructure.leftComp(rxnToRemove,:)=[];
        balanceStructure.rightComp(rxnToRemove,:)=[];
    end
end

newModel=model;
end
