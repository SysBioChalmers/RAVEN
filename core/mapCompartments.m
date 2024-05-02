function geneScoreStructure=mapCompartments(geneScoreStructure,varargin)
% mapCompartments
%   Maps compartments in the geneScoreStructure. This is used if you do not
%   want a models that uses all of the compartment from the predictor. This
%   function will then let you define rules on how the compartments should
%   be merged.
%
%   Any number of rules could be defined as consecutive strings or in a cell array.
%   'comp1'             comp1 should be kept in the structure

%   'comp1=comp2'       The scores in comp2 are merged to comp1 and comp2 is
%                       removed from the structure. This automatically
%                       keeps comp1 in the structure

%   'comp1=comp2 comp3' The scores in comp2 and comp3 are merged to comp1
%                       and comp2 & comp2 are removed from the structure.
%                       This automatically keeps comp1 in the structure

%   'comp1 comp2=comp3' The scores in comp3 are split between comp1 and
%                       comp2. This automatically keeps comp1 and comp2 in
%                       the structure

%   'comp1=other'       The scores in any compartment not included are
%                       merged to comp1. This is applied after all other
%                       rules.
%
%   When one compartment is merged to another the resulting scores will be
%   the best for each gene in either of the compartments. In the case where
%   one compartment is split among several, the scores for the compartment
%   to be merged is weighted with the number of compartments to split to.
%
%   Example: The predictor you use give prediction for Extracellular,
%   Cytosol, Nucleus, Peroxisome, Mitochondria, ER, and Lysosome. You want to
%   have a model with Extracellular, Cytosol, Mitochondria, and Peroxisome
%   where Lysosome is merged with Peroxisome and all other compartments
%   are merged to the Cytosol.
%
%   GSS=mapCompartments(GSS,'Extracellular','Mitochondria','Peroxisome=Lyso
%   some','Cytosol=other');
%
%   geneScoreStructure  a structure to be used in predictLocalization
%
%   Usage: geneScoreStructure=mapCompartments(geneScoreStructure,varargin)

varargin=upper(varargin);

%First find the compartment that will end up in the final structure. They
%are the ones that stand alone or are to the left of some '='
toKeep={};
toMerge={};
I=regexp(varargin,'=','split');
for i=1:numel(varargin)
    if numel(I{i})==1
        toKeep=[toKeep;I{i}];
    else
        J=regexp(I{i}(1),' ','split');
        K=regexp(I{i}(2),' ','split');
        toKeep=[toKeep;J{1}(:)];
        toMerge=[toMerge;K{1}(:)];
    end
end

%Check that there are no compartment that should both be merged and kept
if ~isempty(intersect(toKeep,toMerge))
    EM='There are inconsistencies where one or more compartment(s) should be both kept and merged to another';
    dispEM(EM);
end

%Check that there are no compartments in the rules that are not in the
%geneScoreStructure.
uComps=upper(geneScoreStructure.compartments);
J=[uComps,{'OTHER'}];

if ~isempty(setdiff([toKeep;toMerge],J))
    EM='There are compartment in the rules that are not in geneScoreStructure.compartments';
    dispEM(EM);
end

%Loop through it again and do the mapping
otherIndex=[]; %This stores the rule which maps 'other'.

for i=1:numel(I)
    if numel(I{i})>1
        %Get the compartment indexes that should be mapped
        J=regexp(I{i}(2),' ','split');
        if strcmpi(J{1},'other')
            otherIndex=i;
            continue;
        end
        [k, K]=ismember(J{1},uComps);
        
        %And to where they should be mapped
        J=regexp(I{i}(1),' ','split');
        [l, L]=ismember(J{1},uComps);
        
        %It's not allowed to have rules like A B=C D
        if numel(K)>1 && numel(L)>1
            EM='It is not allowed to have rules like "A B=C D" (map more than one compartment to more than one compartment)';
            dispEM(EM);
        end
        
        if ~all(k) || ~all(l)
            EM='Error in mapping. This most likely means that some compartment(s) are mapped to different compartments in different rules. Use A B=C if you want to map C to several compartments';
            dispEM(EM);
        end
        
        %Get the sum of the scores for the compartments that should be
        %merged to something else
        S=max(geneScoreStructure.scores(:,K),[],2);
        for j=1:numel(L)
            %If the scores are mapped to several different compartments
            %then split the scores between them
            geneScoreStructure.scores(:,L(j))=max(geneScoreStructure.scores(:,L(j)),S./numel(L));
        end
        
        %Remove the comparement that were merged
        geneScoreStructure.compartments(K)=[];
        geneScoreStructure.scores(:,K)=[];
        uComps(K)=[];
    end
end

%Then check if there are remaining compartments that should be removed or
%mapped as 'other'
J=find(~ismember(uComps,toKeep));
if any(J)
    if any(otherIndex)
        K=regexp(I{otherIndex}(1),' ','split');
        [l, L]=ismember(K{1},uComps);
        if l==1 && numel(l)==1
            S=max(geneScoreStructure.scores(:,J),[],2);
            geneScoreStructure.scores(:,L)=max(geneScoreStructure.scores(:,L),S);
        else
            EM='Could not map "other" to more than one compartment';
            dispEM(EM);
        end
    else
        EM='There are compartments that are not defined if they should be kept or removed. Use "A=other" or define more rules if you do not want them to be deleted';
        dispEM(EM,false);
    end
    
    %Remove the comparement that were merged
    geneScoreStructure.compartments(J)=[];
    geneScoreStructure.scores(:,J)=[];
end

%Renormalize
I=max(geneScoreStructure.scores,[],2);
geneScoreStructure.scores=bsxfun(@times, geneScoreStructure.scores, 1./I);

%If there are genes that have score 0 in all compartments, remove them and
%print a warning.
I=find(isnan(geneScoreStructure.scores(:,1))); %Only looks a the first colum as it will be the same for the other ones
if any(I)
    EM='The following genes had score 0.0 in all compartments. They have been removed from the structure. Consider using more rules or "A=other" in order to prevent this:';
    dispEM(EM,false,geneScoreStructure.genes(I));
    geneScoreStructure.scores(I,:)=[];
    geneScoreStructure.genes(I)=[];
end
end
