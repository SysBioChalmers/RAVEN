function [reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,...
    deleteUnconstrained, deleteDuplicates, deleteZeroInterval, deleteInaccessible, deleteMinMax, groupLinear, reservedRxns, suppressWarnings)
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
%   groupLinear           group linear pathways (opt, default false)
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
%           deleteInaccessible, deleteMinMax, groupLinear, reservedRxns,...
%           suppressWarnings)
%
%   Rasmus Agren, 2013-08-01
%

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
    reservedRxns=[];
end
if nargin<9
    suppressWarnings=false;
end

reducedModel=model;
deletedReactions={};
deletedMetabolites={};

if deleteUnconstrained==true
    if isfield(reducedModel,'unconstrained')
        %Remove unbalanced metabolites
        deletedMetabolites=reducedModel.mets(reducedModel.unconstrained~=0);
        reducedModel=removeMets(reducedModel,reducedModel.unconstrained~=0);
        reducedModel=rmfield(reducedModel,'unconstrained');
    end
end

if deleteDuplicates==true
    %Delete all but the last occurrence of duplicate reactions. The
    %reactions must have the same bounds, reversibility, and objective
    %coefficient to be regarded as duplicate
    [reducedModel rxnsToDelete]=contractModel(reducedModel);
    deletedReactions=[deletedReactions; rxnsToDelete];
end

if deleteZeroInterval==true
    %Find all reactions where both LB and UB are 0
    zeroIntervalReactions=and(reducedModel.lb==0,reducedModel.ub==0);
    
    rxnsToDelete=setdiff(reducedModel.rxns(zeroIntervalReactions),reservedRxns);
    deletedReactions=[deletedReactions; rxnsToDelete];   
    
    %Remove reactions
    reducedModel=removeRxns(reducedModel,rxnsToDelete);
    
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
       dispEM('Removing dead-end reactions before removing exchange metabolites',false); 
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
        
        %Also remove metabolites that only participate in one reversible reaction
        %Don't remove the ones in one reversible reaction if is also has a
        %non-zero coefficient in model.b
        notInUse=onlyProducts | onlyReactants | (sum(abs(reducedModel.S')>0)<=1 & (~in & ~out)');
        deletedRxn=false;
        if any(notInUse)        
            %Find their corresponding reactions
            rxnsNotInUse=sum(abs(reducedModel.S(notInUse,:))>0,1)>0;
            rxnsToDelete=setdiff(reducedModel.rxns(rxnsNotInUse),reservedRxns);
            deletedReactions=[deletedReactions; rxnsToDelete];
            
            %Remove reactions
            reducedModel=removeRxns(reducedModel,rxnsToDelete);

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
    reducedModel=removeRxns(reducedModel,rxnsToDelete);
            
    %Remove metabolites
    notInUse=sum(reducedModel.S~=0,2)==0;
    deletedMetabolites=[deletedMetabolites; reducedModel.mets(notInUse)];
    reducedModel=removeMets(reducedModel,notInUse);
end

%Checks that all reactions are irreversible. Might be fixed in the future.
if groupLinear==true
    if ~any(reducedModel.rev)
        if  suppressWarnings==false
            fprintf('NOTE: You have chosen to group linear reactions. This option does not keep track of gene/reaction associations when reactions are merged. Deleting all gene information\n');
        end
        bannedIndexes=getIndexes(reducedModel,reservedRxns,'rxns');
        while 1
        	%Select all metabolites that are only present as reactants/products
        	%in one reaction
            singleNegative=find(sum(reducedModel.S'<0)==1);
            singlePositive=find(sum(reducedModel.S'>0)==1);
   
            %Retrieve the common metabolites
            common=intersect(singleNegative,singlePositive);
            
            %Merge the reactions for each of these metabolites
            mergedSomeRxn=false;
            for i=1:numel(common)
            	involvedRxns=find(reducedModel.S(common(i),:));
                
                %Check so that it doesn't try to merge reactions that are
                %restricted. The second check is needed because it could be
                %that some of the metabolites have already been merged in
                %this iteration so that the reactions no longer exist
                if isempty(intersect(bannedIndexes,involvedRxns)) && any(involvedRxns)
                	%This is the ratio that described how many times the second
                    %reaction should be multiplied before being merged with
                    %the first
                    stoichRatio=abs(reducedModel.S(common(i),involvedRxns(1))/reducedModel.S(common(i),involvedRxns(2)));
                    reducedModel.S(:,involvedRxns(1))=reducedModel.S(:,involvedRxns(1))+reducedModel.S(:,involvedRxns(2))*stoichRatio;
                    reducedModel.S(common(i),involvedRxns(1))=0;
                    
                    %Change the name and id to reflect the merging
                    reducedModel.rxns{involvedRxns(1)}=[reducedModel.rxns{involvedRxns(1)} '_' reducedModel.rxns{involvedRxns(2)}];
                    reducedModel.rxnNames{involvedRxns(1)}=reducedModel.rxns{involvedRxns(1)};
                    
                    %Remove the second reaction by setting all coefficients
                    %to 0
                    reducedModel.S(:,involvedRxns(2))=0;
                    mergedSomeRxn=true;
                end
            end
            
            %If it couldn't merge anything then exit the loop
            if mergedSomeRxn==false
                break;
            end
        end
        
        %Now delete all reactions that involve no metabolites
        I=sum(reducedModel.S~=0)==0;
        
        %Remove reactions
        deletedReactions=[deletedReactions; reducedModel.rxns(I)];
        reducedModel=removeRxns(reducedModel,find(I));
            
        %Remove metabolites
        notInUse=sum(reducedModel.S~=0,2)==0;
        deletedMetabolites=[deletedMetabolites; reducedModel.mets(notInUse)];
        reducedModel=removeMets(reducedModel,notInUse);
        
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
    else
    	dispEM('Cannot group reactions for a model that has reversible reactions');
    end
end

end
