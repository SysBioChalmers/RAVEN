function [model, deletedRxns, duplicateRxns]=mergeCompartments(model,varargin)
% mergeCompartments  Merge all compartments in a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% keepUnconstrained : logical, optional
%     keep metabolites that are unconstrained in a 'unconstrained'
%     compartment. If these are merged the exchange reactions will most often
%     be deleted (default false).
% deleteRxnsWithOneMet : logical, optional
%     delete reactions with only one metabolite. These reactions come from
%     reactions such as A[c] + B[c] => A[m]. In some models hydrogen is
%     balanced around each membrane with reactions like this (default
%     false).
% distReverse : logical, optional
%     distinguish reactions with same metabolites but different reversibility
%     as different reactions (default true).
%
% Returns
% -------
% model : struct
%     a model with all reactions located to one compartment.
% deletedRxns : cell
%     reactions that were deleted because of only having one metabolite
%     after merging.
% duplicateRxns : cell
%     identical reactions that occurred in different compartments and were
%     deleted because they turned to be duplicated after merging.
%
% Examples
% --------
%     [model, deletedRxns, duplicateRxns] = mergeCompartments(model);
%
% Notes
% -----
% Merges all compartments into one 's' compartment (for 'System'). This can
% be useful for example to ensure that there are metabolic capabilities to
% synthesize all metabolites.
%
% If the metabolite IDs reflect the compartment that they are in the IDs may
% no longer be representative.

p=parseRAVENargs(varargin, {'keepUnconstrained',false; 'deleteRxnsWithOneMet',false; 'distReverse',true});
keepUnconstrained=p.keepUnconstrained;
deleteRxnsWithOneMet=p.deleteRxnsWithOneMet;
distReverse=p.distReverse;

if ~isfield(model,'unconstrained')
    keepUnconstrained=false;
end

%Keep track of which reactions only contained one metabolite before the
%merging as they are probable exchange reactions.
if deleteRxnsWithOneMet==true
    reservedRxns=model.rxns(sum(model.S~=0)==1);
    if ~isempty(reservedRxns) && isfield(model,'unconstrained')
        %If there is no unconstrained field these reactions are probably
        %exchange reactions and shall be kept. If not then print a warning
        EM='There are reactions with only one metabolite. Cannot determine whether they are exchange reactions since there is no unconstrained field';
        dispEM(EM,false);
    end
end

%Loop through each metabolite, and if it is not unconstrained then change
%the S matrix to use the metabolite with the lowest index in model.comps
%instead
uNames=unique(model.metNames);
for i=1:numel(uNames)
    %Find all metabolites with this name..
    I=ismember(model.metNames,uNames(i));
    
    %Find the first of those that is not unconstrained. This is the one
    %that the other "un-unconstrained" should be changed to.
    if keepUnconstrained==true
        mergeTo=find(I & model.unconstrained==false,1);
        
        %This could happen if there is only one metabolite and it is
        %unconstrained
        if isempty(mergeTo)
            continue;
        end
    else
        mergeTo=find(I,1);
    end
    I(mergeTo)=false; %Do not do anything for itself
    I=find(I);
    
    %Go through each of the metabolites with this name and update them to
    %mergeTo
    for j=1:numel(I)
        if keepUnconstrained==true && model.unconstrained(I(j))==true
            continue;
        end
        %Update the S matrix
        model.S(mergeTo,:)=model.S(mergeTo,:)+model.S(I(j),:);
        model.S(I(j),:)=0;
    end
end

%Remove all metabolites that are no longer used (with a bit of a trick)
model=removeReactions(model,{},true);

%Update all mets so that they are in compartment "s" with id "1"
model.compNames={'System'};
model.comps={'s'};
model.compOutside={''};
model.metComps=ones(numel(model.mets),1);

%Add exchange mets to another compartment "b" with id "2"
if keepUnconstrained==true
    model.compNames=[model.compNames; {'Unconstrained'}];
    model.comps=[model.comps; {'b'}];
    model.compOutside=[model.compOutside; {'s'}];
    model.metComps(model.unconstrained~=0)=2;
end

%Transport reactions on the form A <=> B will have been deleted by the
%merging. Remove those reactions
model=removeMets(model,{},false,true,true);

if deleteRxnsWithOneMet==true
    I=model.rxns(sum(model.S~=0)==1);
    
    %Delete the reactions that contain only one metabolite after the
    %merging but not before
    deletedRxns=setdiff(I,reservedRxns);
    model=removeReactions(model,deletedRxns,true,true);
else
    deletedRxns={};
end

%And then finally merge the identical reactions
[model, ~, duplicateRxns]=contractModel(model,distReverse);
end
