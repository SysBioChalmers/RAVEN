function [notProduced, notProducedNames, neededForProductionMat,minToConnect,model]=checkProduction(model,checkNeededForProduction,excretionFromCompartments,printDetails)
% checkProduction
%   Checks which metabolites that can be produced from a model using the
%   specified constraints.
%
%   model                       a model structure
%   checkNeededForProduction    for each of the metabolites that could not
%                               be produced, include an artificial
%                               production reaction and calculate which new
%                               metabolites that could be produced as en
%                               effect of this (opt, default false)
%   excretionFromCompartments   cell array with compartment ids from which
%                               metabolites can be excreted (opt, default
%                               model.comps)
%   printDetails                print details to the screen (opt, default
%                               true)
%
%   notProduced                 cell array with metabolites that could not
%                               be produced
%   notProducedNames            cell array with names and compartments for
%                               metabolites that could not be produced
%   neededForProductionMat      matrix where n x m is true if metabolite n
%                               allows for production of metabolite m
%   minToConnect                structure with the minimal number of
%                               metabolites that need to be connected in
%                               order to be able to produce all other
%                               metabolites and which metabolites each of
%                               them connects
%   model                       updated model structure with excretion
%                               reactions added
%
%   The function is intended to be used to identify which metabolites must
%   be connected in order to have a fully connected network. It does so by
%   first identifying which metabolites could have a net production in the
%   network. Then it calculates which other metabolites must be able to
%   have net production in order to have production of all metabolites in
%   the network. So, if a network contains the equations A[external]->B,
%   C->D, and D->E it will identify that production of C will connect
%   the metabolites D and E.
%
%   Usage: [notProduced, notProducedNames,neededForProductionMat,minToConnect,model]=...
%           checkProduction(model,checkNeededForProduction,...
%           excretionFromCompartments,printDetails)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<2
    checkNeededForProduction=false;
end

if nargin<3
    excretionFromCompartments=model.comps;
end

if nargin<4
    printDetails=true;
end

%Add an exchange reaction for each metabolite in the allowed compartments
%and see if it can carry a flux
allowedMetIds=ismember(model.comps(model.metComps),excretionFromCompartments);
allowedMetIndexes=find(allowedMetIds);
[model, addedRxns]=addExchangeRxns(model,'out',allowedMetIndexes);

canProduce=haveFlux(model,10^-5,addedRxns);

notProduced=find(~canProduce);
minToConnect={};
if checkNeededForProduction==true
    %For each of the metabolites that couldn't be produced allow uptake and
    %check which of the other metabolites that couldn't be produced that
    %can be produced
    neededForProductionMat=false(numel(notProduced));
    for i=1:numel(notProduced)
        %Add uptake for this metabolite
        if i>1
            %Reset last iteration
            model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i-1))=model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i-1))*-1;
        end
        %Change the production reaction to an uptake reaction
        model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i))=model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i))*-1;
        
        %Test which of the metabolites that couldn't be produced that can
        %be produced now
        neededForProductionMat(i,:)=haveFlux(model,10^-5,addedRxns(notProduced));
    end
    %Calculate the smallest number of metabolites that must be connected to
    %make everything connected and return their names
    
    %The algorithm is relatively straight forward. It finds the metabolite
    %that connects the most unconnected metabolites (iteratively), adds it
    %and removes the now connected metabolites until all are connected.
    %This is not guaranteed to find the global minimum
    neededForProdTemp=neededForProductionMat;
    while 1==1
        %Get all metabolites a metabolite makes connected
        totalConnected=false(size(neededForProdTemp));
        for i=1:numel(notProduced)
            totalConnected(i,:)=neededForProdTemp(i,:);
            
            lastIter=0;
            while 1==1
                [~, a]=find(neededForProdTemp(totalConnected(i,:),:));
                totalConnected(i,a)=true;
                if numel(a)==lastIter
                    break; %No more connections were possible
                else
                    lastIter=numel(a);
                end
            end
        end
        [connections, mostConnected]=max(sum(totalConnected,2));
        
        if connections>0
            %Add the most connected metabolite to the list and remove all
            %metabolites that it's connected to
            metID=allowedMetIndexes(notProduced(mostConnected));
            entry=[model.metNames{metID},'[',model.comps{model.metComps(metID)},'] (connects ' num2str(connections) ' metabolites)'];
            minToConnect=[minToConnect;entry];
            neededForProdTemp(totalConnected(mostConnected,:),:)=false;
            neededForProdTemp(:,totalConnected(mostConnected,:))=false;
        else
            break;
        end
    end
else
    neededForProductionMat=[];
end

notProducedNames=strcat(model.metNames(allowedMetIndexes(notProduced)),'[',model.comps(model.metComps(allowedMetIndexes(notProduced))),']');

if printDetails==true
    fprintf('The following metabolites could not be produced:\n');
    [notProducedNamesTemp,perm]=sort(notProducedNames);
    
    if checkNeededForProduction==true
        neededForProdTemp=neededForProductionMat(:,perm);
        neededForProdTemp=neededForProdTemp(perm,:);
        fprintf('\tIf the production of a metabolite is dependent on some other metabolites then those are printed under the name\n\n');
    end
    for i=1:numel(notProducedNamesTemp)
        fprintf([notProducedNamesTemp{i} '\n']);
        neededForProdTemp(i,i)=false; %Not neat to do this here. Prevent printing itself
        if checkNeededForProduction==true
            enablesProduction=find(neededForProdTemp(:,i));
            if any(enablesProduction)
                for j=1:numel(enablesProduction)
                    fprintf(['\t' notProducedNamesTemp{enablesProduction(j)} '\n']);
                end
            end
        end
    end
end
end
