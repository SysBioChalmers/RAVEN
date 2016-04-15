function [model deletedRxns]=mergeCompartments(model,keepUnconstrained,deleteRxnsWithOneMet)
% mergeCompartments
%   Merge all compartments in a model
%
%   model                 a model structure
%   keepUnconstrained     keep metabolites that are unconstrained in a
%                         'unconstrained' compartment. If these are merged the
%                         exchange reactions will most often be deleted (opt,
%                         default false)
%   deleteRxnsWithOneMet  delete reactions with only one metabolite. These
%                         reactions come from reactions such as A[c] + B[c]
%                         => A[m]. In some models hydrogen is balanced around
%                         each membrane with reactions like this (opt,
%                         default false)
%
%   model                 a model with all reactions located to one compartment
%   deletedRxns           reactions that were deleted because of only
%                         having one metabolite after merging 
%
%   Merges all compartments into one 's' compartment (for 'System'). This can
%   be useful for example to ensure that there are metabolic capabilities to
%   synthesize all metabolites.
%
%   NOTE: If the metabolite IDs reflect the compartment that they are in
%   the IDs may no longer be representative.
%
%   Usage: model=mergeCompartments(model,keepUnconstrained,deleteRxnsWithOneMet)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    keepUnconstrained=false;
end
if nargin<3
    deleteRxnsWithOneMet=false;
end

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
        fprintf('WARNING: There are reactions with only one metabolite. Cannot determine whether they are exchange reactions since there is no unconstrained field.\n');
    end
end

%Loop through each metabolite, and if it's not unconstrained then change
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
        
        %This could happen if there is only one metabolite and it's
        %unconstrained
        if isempty(mergeTo)
            continue;
        end
    else
        mergeTo=find(I,1);
    end
    I(mergeTo)=false; %Don't do anything for itself
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
model=removeRxns(model,{},true);

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
    model=removeRxns(model,deletedRxns,true,true);
else
    deletedRxns={};
end

%And then finally merge the identical reactions
model=contractModel(model);
end
