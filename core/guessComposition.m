function [model, guessedFor, couldNotGuess]=guessComposition(model, printResults)
% guessComposition
%   Attempts to guess the composition of metabolites without information
%   about elemental composition
%
%   model               a model structure
%   printResults        true if the output should be printed (opt, default true)
%
%   model               a model structure with information about elemental
%                       composition added
%   guessedFor          indexes for the metabolites for which a composition
%                       could be guessed
%   couldNotGuess       indexes for the metabolites for which no
%                       composition could be assigned
%
%   This function works in a rather straight forward manner:
%
%   1. Get the metabolites which lack composition and participates in
%   at least one reaction where all other metabolites have composition information
%   2. Loop through them and calculate their composition based on the rest
%   of the involved metabolites. If there are any inconsistencies, so that
%   a given metabolite should have different composition in different
%   equations, then throw an error
%   3. Go to 1
%
%   This simple approach requires that the rest of the metabolites have
%   correct composition information, and that the involved reactions are
%   correct. The function will exit with an error on any inconsistencies,
%   which means that it could also be used as a way of checking the model
%   for errors. Note that just because this exits sucessfully, the
%   calculated compositions could still be wrong (in case that the existing
%   compositions were wrong)
%
%   Usage: [newModel, guessedFor, couldNotGuess]=guessComposition(model, printResults)
%
%   Rasmus Agren, 2013-11-06
%

if nargin<2
    printResults=true;
end

%The metabolites for which there is no elemental composition
originalMissing=unique(model.metNames(cellfun(@isempty,model.metFormulas)));
guessedFor={};

%This is to keep track of if new metabolite compositions were predicted
predicted=true;
while predicted==true
    predicted=false;

    %Get the unique names (composition should be independent of compartment)
    %for the metabolites without composition
    metNames=unique(model.metNames(cellfun(@isempty,model.metFormulas)));

    %Parse the formulas in the model
    [elements, useMat, exitFlag]=parseFormulas(model.metFormulas, true,false);

    for i=1:numel(metNames)
        %Get the metabolites with this name. Not so neat, but this is a
        %fast function anyways
        mets=find(ismember(model.metNames,metNames(i)));
        
        currentComp=[];

        %Loop through the metabolites
        %-1: Could not assign due to missing info. -2: Could not assign due to contradiction
        %1: Composition assigned
        metStatus=-1;
        for j=1:numel(mets)
            %Get the reactions that the metabolite participates in
            [crap, I]=find(model.S(mets(j),:));
            if any(I)
                for k=1:numel(I)
                    %Loop through the reactions and check if all other mets in them
                    %have known composition
                    eqn=model.S(:,I(k));
                    eqn(mets(j))=0;
                    if all(exitFlag(eqn~=0)==1)
                        %This means that all other mets had composition. Calculate
                        %the resulting composition for the unknown one
                        comp=useMat'*eqn;
                        
                        %This can result in round off errors if there are
                        %stoichiometries with many decimals. Ignore values
                        %below 10^-12
                        comp(abs(comp)<10^-12)=0;

                        %Check if the composition consist of both negative and
                        %positive values. If so, throw an error
                        if all(comp<=0) || all(comp>=0)
                            comp=abs(comp);
                            if isempty(currentComp)
                                currentComp=comp;
                            end
                            %Also to deal with round off errors
                            if all(abs(currentComp-comp)<10^-10)
                                metStatus=1;
                            else
                                metStatus=-2;
                                break;

                                %%Check if there is an inconcistency
                                %if any(currentComp~=comp)
                                %    dispEM(['Could not predict composition of ' model.metNames{mets(i)} ],false);
                                %end
                            end
                        else
                            %Check if there is an inconcistency
                            %if any(currentComp~=comp)
                            %    dispEM(['Could not predict composition of ' model.metNames{loopThrough(i)} ],false);
                            %end
                            metStatus=-2;
                            break;
                        end
                    end
                end
                %If there was contradictions in one compartment, then abort for
                %all compartments
                if metStatus==-2
                    break;
                end
            else
                %The metabolite doesn't participate, no composition can be
                %calculated
                metStatus=-1;
            end
        end
        %Check status of the metabolite
        switch metStatus
            case -2
                dispEM(['Could not predict composition for "' metNames{i} '" due to inconsistencies'],false);
            case 1
                %Calculate and add the composition
                str=getCompString(elements,comp);
                model.metFormulas(mets)={str};
                if printResults==true
                    fprintf(['Predicted composition for "' metNames{i} '" to be ' str '\n']);
                end
                
                %Keep track
                guessedFor=[guessedFor;metNames(i)];
                
                predicted=true; %To loop again
        end        
    end
end

couldNotGuess=setdiff(originalMissing,guessedFor);
end

%Helper function for getting the composition string
function str=getCompString(elements,comp)
    str='';
    
    for i=1:numel(comp)
       if comp(i)~=0
          if comp(i)==1
             str=[str  elements.abbrevs{i}];
          else
             str=[str  elements.abbrevs{i} num2str(comp(i))]; 
          end
       end
    end
end