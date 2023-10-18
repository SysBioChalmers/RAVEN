function model=sortModel(model,sortReversible,sortMetName,sortReactionOrder)
% sortModel
%   Sorts a model based on metabolite names and compartments
%
%   model             a model structure
%   sortReversible    sorts the reversible reactions so the the metabolite
%                     that is first in lexiographical order is a reactant
%                     (opt, default true)
%   sortMetName       sort the metabolite names in the equation, also uses
%                     compartment abbreviation (opt, default false)
%   sortReactionOrder sorts the reaction order within each subsystem so that
%                     reactions consuming some metabolite comes efter
%                     reactions producing it. This overrides the
%                     sortReversible option and reactions are sorted so that
%                     the production direction matches the consumption
%                     direction (opt, default false)
%
%   model             an updated model structure
%
%   Usage: model=sortModel(model,sortReversible,sortMetName,sortReactionOrder)

if nargin<2
    sortReversible=true;
end
if nargin<3
    sortMetName=false;
end
if nargin<4
    sortReactionOrder=false;
end

if sortMetName==true
    %Assuming that metComps are the indexes. Should be changed at one point
    [~, metIndexes]=sort(strcat(model.metNames,'[',model.comps(model.metComps),']'));
    model=permuteModel(model,metIndexes,'mets');
end

if sortReversible==true && sortReactionOrder==false
    %Get all reversible reactions
    revIndexes=find(model.rev);
    
    %Loop through them
    for i=1:numel(revIndexes)
        %Create a cell array with all the metabolite names
        mets=find(model.S(:,revIndexes(i)));
        metNames=strcat(model.metNames(mets),model.comps(model.metComps(mets)));
        
        if iscellstr(metNames)
            [~, indexes]=sort(metNames);
            
            if model.S(mets(indexes(1)),revIndexes(i))>0
                model.S(:,revIndexes(i))=model.S(:,revIndexes(i))*-1;
            end
        end
    end
end

if sortReactionOrder==true
    %Check if the model has sub-systems, otherwise throw an error
    if ~isfield(model,'subSystems')
        EM='The model must contain a subSystems field in order to sort reaction order';
        dispEM(EM);
    end
    
    subsystemsUnique='';
    subsystemsConcatenated='';
    for i=1:numel(model.subSystems)
        subsystemsConcatenated{i,1}=strjoin(model.subSystems{i,1},';');
        if ~isempty(model.subSystems{i,1})
            for j=1:numel(model.subSystems{i,1})
                subsystemsUnique{numel(subsystemsUnique)+1,1}=model.subSystems{i,1}{1,j};
            end
        end
    end
    subsystemsUnique=unique(subsystemsUnique);
    for i=1:numel(subsystemsUnique)
        %Get all reactions for that subsystem
        rxns=find(~cellfun(@isempty,regexp(subsystemsConcatenated,subsystemsUnique(i))));
        
        %Temporarily ignore large subsystems because of inefficient
        %implementation
        if numel(rxns)<2 || numel(rxns)>250
            continue;
        end
        
        nRxns=numel(rxns);
        revRxns=rxns(model.rev(rxns)~=0);
        
        %This variable will hold the current reversibility directions of
        %the reversible reactions. 1 means the same direction as in the
        %original model and -1 means the opposite direction.
        oldRev=ones(numel(revRxns),1);
        
        %The problem could probably be solved analytically but a simple
        %random method is implemented here instead. Two reactions are
        %chosen randomly and their positions are switched. A score is
        %calculated based on the number of metabolites that are produced
        %before they are consumed. If the perturbed model has a better or
        %equal score than the original the reaction order is switched. If
        %no increase in score has been seen after 1000*rxnsInSubsystem then
        %the optimization is terminated
        
        rxnOrder=1:nRxns;
        oldScore=-inf;
        counter=0;
        firstIter=true;
        while 1==1
            counter=counter+1;
            if counter==100*nRxns
                break;
            end
            
            newRxnOrder=rxnOrder;
            rev=oldRev;
            
            if firstIter==false
                y=randperm(nRxns,2);
                
                %Switch the order
                newRxnOrder(y(1))=rxnOrder(y(2));
                newRxnOrder(y(2))=rxnOrder(y(1));
                
                %With a 50% chance, also switch the reversibility of one of
                %the reactions
                if rand()>0.5 && numel(rev)>1
                    n=randperm(numel(rev),1);
                    rev(n)=rev(n)*-1;
                end
            end
            firstIter=false;
            
            tempS=model.S;
            
            %Switch the directionalities
            for j=1:numel(rev)
                if rev(j)==-1
                    tempS(:,revRxns(j))=tempS(:,revRxns(j)).*-1;
                end
            end
            
            %Get the metabolites that are involved and when they are
            %produced/consumed
            s=tempS(:,newRxnOrder);
            
            %Remove mets that aren't used in both directions
            s=s(any(s,2),:);
            
            %Add so that all mets are produced and consumed in the end
            s=[s ones(size(s,1),1) ones(size(s,1),1)*-1];
            
            %For each metabolite, find the reaction where it's first
            %produced and the reaction where it's first consumed
            s1=s>0;
            r1=arrayfun(@(x) find(s1(x,:),1,'first'),1:size(s1,1));
            s2=s<0;
            r2=arrayfun(@(x) find(s2(x,:),1,'first'),1:size(s2,1));
            
            score=sum(r1<r2);
            
            if score>=oldScore
                if score>oldScore
                    counter=0;
                end
                oldScore=score;
                oldRev=rev;
                rxnOrder=newRxnOrder;
            end
        end
        
        %Update the model for this subsystem
        for j=1:numel(oldRev)
            if oldRev(j)==-1
                model.S(:,revRxns(j))=model.S(:,revRxns(j)).*-1;
            end
        end
        order=1:numel(model.rxns);
        order(rxns)=rxns(rxnOrder);
        model=permuteModel(model, order, 'rxns');
    end
end
end
