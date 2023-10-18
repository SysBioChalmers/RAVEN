function [reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,...
    deleteUnconstrained, deleteDuplicates, deleteZeroInterval, deleteInaccessible, deleteMinMax, groupLinear, constrainReversible, reservedRxns, suppressWarnings)
% simplifyModel
%   Simplifies a model by deleting reactions/metabolites
%
%   model                 a model structure
%   deleteUnconstrained   delete metabolites marked as unconstrained (opt, default true)
%   deleteDuplicates      delete all but one of duplicate reactions (opt, default false)
%   deleteZeroInterval    delete reactions that are constrained to zero flux (opt, default false)
%   deleteInaccessible    delete dead end reactions (opt, default false)
%   deleteMinMax          delete reactions that cannot carry a flux by trying
%                         to minimize/maximize the flux through that
%                         reaction. May be time consuming (opt, default false)
%   groupLinear           group linearly dependent pathways (opt, default false)
%   constrainReversible   check if there are reversible reactions which can
%                         only carry flux in one direction, and if so
%                         constrain them to be irreversible. This tends to
%                         allow for more reactions grouped when using
%                         groupLinear (opt, default false)
%   reservedRxns          cell array with reaction IDs that are not allowed to be
%                         removed (opt)
%   suppressWarnings      true if warnings should be suppressed (opt,
%                         default false)
%
%   reducedModel          an updated model structure
%   deletedReactions      a cell array with the IDs of all deleted reactions
%   deletedMetabolites    a cell array with the IDs of all deleted
%                         metabolites
%
%   This function is for reducing the model size by removing
%   reactions and associated metabolites that cannot carry flux. It can also
%   be used for identifying different types of gaps.
%
%   Usage: [reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,...
%           deleteUnconstrained, deleteDuplicates, deleteZeroInterval,...
%           deleteInaccessible, deleteMinMax, groupLinear,...
%           constrainReversible, reservedRxns, suppressWarnings)

if nargin<2
    deleteUnconstrained=true;
end
if nargin<3
    deleteDuplicates=false;
end
if nargin<4
    deleteZeroInterval=false;
end
if nargin<5
    deleteInaccessible=false;
end
if nargin<6
    deleteMinMax=false;
end
if nargin<7
    groupLinear=false;
end
if nargin<8
    constrainReversible=false;
end
if nargin<9
    reservedRxns=[];
else
    reservedRxns=convertCharArray(reservedRxns);
end
if nargin<10
    suppressWarnings=false;
end

reducedModel=model;
deletedReactions={};
deletedMetabolites={};

if deleteUnconstrained==true
    if isfield(reducedModel,'unconstrained')
        %Remove unbalanced metabolites
        deletedMetabolites=reducedModel.mets(reducedModel.unconstrained~=0);
        reducedModel=removeMets(reducedModel,reducedModel.unconstrained~=0,false,false,false,true);
        reducedModel=rmfield(reducedModel,'unconstrained');
    end
end

if deleteDuplicates==true
    %Delete all but the last occurrence of duplicate reactions. The
    %reactions must have the same bounds, reversibility, and objective
    %coefficient to be regarded as duplicate
    [reducedModel, rxnsToDelete, ~]=contractModel(reducedModel);
    deletedReactions=[deletedReactions; rxnsToDelete];
end

if deleteZeroInterval==true
    %Find all reactions where both LB and UB are 0
    zeroIntervalReactions=and(reducedModel.lb==0,reducedModel.ub==0);
    
    rxnsToDelete=setdiff(reducedModel.rxns(zeroIntervalReactions),reservedRxns);
    deletedReactions=[deletedReactions; rxnsToDelete];
    
    %Remove reactions
    reducedModel=removeReactions(reducedModel,rxnsToDelete);
    
    %Find metabolites that no longer are used and delete them
    notInUse=sum(reducedModel.S~=0,2)==0;
    deletedMetabolites=[deletedMetabolites;reducedModel.mets(notInUse)];
    
    %Remove metabolites
    reducedModel=removeMets(reducedModel,notInUse);
end

if deleteInaccessible==true
    %Print a warning if exchange metabolites haven't been deleted yet. This
    %often means that the only allowed products are the metabolites that
    %are taken up be the system
    if isfield(reducedModel,'unconstrained') && suppressWarnings==false
        EM='Removing dead-end reactions before removing exchange metabolites';
        dispEM(EM,false);
    end
    
    while true
        %Find all metabolites that are only reactants or only products.
        %Delete those metabolites and reactions. This is done until no more
        %such reactions are found
        
        %Construct a stoichiometric matrix where all reactions are written
        %as reversible
        
        %Add fake exchange reactions if model.b~=0
        in=any(reducedModel.b<0,2);
        out=any(reducedModel.b>0,2);
        I=speye(numel(reducedModel.mets));
        revS=[reducedModel.S,reducedModel.S(:,reducedModel.rev~=0)*-1 I(:,in) I(:,out)*-1];
        
        metUsage=sum(abs(revS')>0);
        onlyProducts=sum(revS'>0) == metUsage;
        onlyReactants=sum(revS'<0) == metUsage;
        
        %Also remove metabolites that only participate in one reversible
        %reaction Don't remove the ones in one reversible reaction if is
        %also has a non-zero coefficient in model.b
        notInUse=onlyProducts | onlyReactants | (sum(abs(reducedModel.S')>0)<=1 & (~in & ~out)');
        deletedRxn=false;
        if any(notInUse)
            %Find their corresponding reactions
            rxnsNotInUse=sum(abs(reducedModel.S(notInUse,:))>0,1)>0;
            rxnsToDelete=setdiff(reducedModel.rxns(rxnsNotInUse),reservedRxns);
            deletedReactions=[deletedReactions; rxnsToDelete];
            
            %Remove reactions
            reducedModel=removeReactions(reducedModel,rxnsToDelete);
            
            %Remove metabolites. Recalculate since it could be that some
            %cannot be deleted due to reserved rxns
            notInUse=sum(reducedModel.S~=0,2)==0;
            deletedMetabolites=[deletedMetabolites; reducedModel.mets(notInUse)];
            reducedModel=removeMets(reducedModel,notInUse);
            
            %It could happen that it just finds reserved reactions and gets
            %stuck
            if ~isempty(rxnsToDelete)
                deletedRxn=true;
            end
        else
            break;
        end
        if deletedRxn==false
            break;
        end
    end
end

if deleteMinMax==true
    %Get reactions that can't carry fluxes. This should be done
    %algebraically if possible
    I=~haveFlux(reducedModel);
    
    %Remove reactions
    rxnsToDelete=setdiff(reducedModel.rxns(I),reservedRxns);
    deletedReactions=[deletedReactions; rxnsToDelete];
    reducedModel=removeReactions(reducedModel,rxnsToDelete);
    
    %Remove metabolites
    notInUse=sum(reducedModel.S~=0,2)==0;
    deletedMetabolites=[deletedMetabolites; reducedModel.mets(notInUse)];
    reducedModel=removeMets(reducedModel,notInUse);
end

if constrainReversible==true
    revs=find(reducedModel.rev);
    [I,J]=getAllowedBounds(reducedModel,revs);
    
    I=abs(I);
    J=abs(J);
    
    %Get the "small" values
    K=I<10^-10;
    L=J<10^-10;
    
    %Keep the values where only one direction is zero (small)
    I(K==L)=[];
    J(K==L)=[];
    revs(K==L)=[];
    
    %Change the reversibility of the remaining reactions
    reducedModel.rev(revs(J>10^-10))=0; %Ignore reverse direction
    reducedModel.lb(revs(J>10^-10))=0;
    
    toSwitch=revs(I>10^-10);
    reducedModel.rev(toSwitch)=0; %Change directionality
    reducedModel.ub(toSwitch)=reducedModel.lb(toSwitch)*-1;
    reducedModel.lb(toSwitch)=0;
    reducedModel.S(:,toSwitch)=reducedModel.S(:,toSwitch).*-1;
    reducedModel.c(toSwitch)=reducedModel.c(toSwitch)*-1;
end
if groupLinear==true
    if  suppressWarnings==false
        fprintf('NOTE: You have chosen to group linear reactions. This option does not keep track of gene/reaction associations when reactions are merged. Deleting all gene information\n');
    end
    
    reducedModel.genes={};
    reducedModel.rxnGeneMat=sparse(numel(reducedModel.rxns),0);
    reducedModel.grRules(:)={''};
    
    if isfield(reducedModel,'geneShortNames')
        reducedModel.geneShortNames={};
    end
    if isfield(reducedModel,'geneMiriams')
        reducedModel.geneMiriams={};
    end
    if isfield(reducedModel,'geneComps')
        reducedModel.geneComps=[];
    end
    
    %Convert the model to irreversible
    irrevModel=convertToIrrev(reducedModel);
    
    %Loop through and iteratively group linear reactions
    while 1
        %Get the banned reaction indexes. Note that the indexes will change
        %in each iteration, but the names will not as they won't be merged
        %with any other reaction
        bannedIndexes=getIndexes(irrevModel,reservedRxns,'rxns');
        
        %Select all metabolites that are only present as reactants/products
        %in one reaction
        singleNegative=find(sum(irrevModel.S'<0)==1);
        singlePositive=find(sum(irrevModel.S'>0)==1);
        
        %Retrieve the common metabolites
        common=intersect(singleNegative,singlePositive);
        
        mergedSome=false;
        
        %Loop through each of them and see if the reactions should be
        %merged
        for i=1:numel(common)
            involvedRxns=find(irrevModel.S(common(i),:));
            
            %Check so that one or both of the reactions haven't been merged
            %already
            if numel(involvedRxns)==2 && isempty(intersect(bannedIndexes,involvedRxns))
                %Calculate how many times the second reaction has to be
                %multiplied before being merged with the first
                stoichRatio=abs(irrevModel.S(common(i),involvedRxns(1))/irrevModel.S(common(i),involvedRxns(2)));
                
                %Add the second to the first
                irrevModel.S(:,involvedRxns(1))=irrevModel.S(:,involvedRxns(1))+irrevModel.S(:,involvedRxns(2))*stoichRatio;
                
                %Clear the second reaction
                irrevModel.S(:,involvedRxns(2))=0;
                
                %This is to prevent numerical issues. It should be 0
                %already
                irrevModel.S(common(i),involvedRxns(1))=0;
                
                %At this point the second reaction is certain to be deleted
                %in a later step and can therefore be ignored
                
                %Recalculate the bounds for the new reaction. This can be
                %problematic since the scale of the bounds may change
                %dramatically. Let the most constraining reaction determine
                %the new bound
                lb1=irrevModel.lb(involvedRxns(1));
                lb2=irrevModel.lb(involvedRxns(2));
                ub1=irrevModel.ub(involvedRxns(1));
                ub2=irrevModel.ub(involvedRxns(2));
                
                if lb2~=-inf
                    irrevModel.lb(involvedRxns(1))=max(lb1,lb2/stoichRatio);
                end
                if ub2~=inf
                    irrevModel.ub(involvedRxns(1))=min(ub1,ub2/stoichRatio);
                end
                
                %Then recalculate the objective coefficient. The resulting
                %coefficient is the weighted sum of the previous
                irrevModel.c(involvedRxns(1))=irrevModel.c(involvedRxns(1))+irrevModel.c(involvedRxns(2))*stoichRatio;
                
                %Iterate again
                mergedSome=true;
            end
        end
        
        %All possible reactions merged
        if mergedSome==false
            break;
        end
        
        %Now delete all reactions that involve no metabolites
        I=find(sum(irrevModel.S~=0)==0);
        
        %Remove reactions
        irrevModel=removeReactions(irrevModel,I);
        
        %Remove metabolites
        notInUse=sum(irrevModel.S~=0,2)==0;
        irrevModel=removeMets(irrevModel,notInUse);
    end
    
    reducedModel=irrevModel;
end
end
