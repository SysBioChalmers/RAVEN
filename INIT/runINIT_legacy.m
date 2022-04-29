function [outModel, deletedRxns, metProduction, fValue]=runINIT_legacy(model,rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,noRevLoops,params)
% runINIT
%	Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks. This is the 
%   original implementation, which is now replaced with ftINIT.
%
%   model           a reference model structure
%   rxnScores       a vector of scores for the reactions in the model.
%                   Positive scores are reactions to keep and negative
%                   scores are reactions to exclude (opt, default all 0.0)
%   presentMets     cell array with unique metabolite names that the model
%                   should produce (opt, default [])
%   essentialRxns   cell array of reactions that are essential and that
%                   have to be in the resulting model. This is normally
%                   used when fitting a model to task (see fitTasks) (opt,
%                   default [])
%   prodWeight      a score that determines the value of having
%                   net-production of metabolites. This is a way of having
%                   a more functional network as it provides a reason for
%                   including bad reactions for connectivity reasons. This
%                   score is for each metabolite, and the sum of these weights
%                   and the scores for the reactions is what is optimized
%                   (opt, default 0.5)
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                   This results in fewer reactions being considered
%                   dead-ends, but all reactions in the resulting model may
%                   not be able to carry flux. If this is "false" then the
%                   equality constraints are taken from model.b. If the
%                   input model lacks exchange reactions then this should
%                   probably be "true", or a large proportion of the model
%                   would be excluded for connectivity reasons
%                   (opt, default false)
%   noRevLoops      true if reversible reactions should be constrained to
%                   only carry flux in one direction. This prevents
%                   reversible reactions from being wrongly assigned as
%                   connected (the forward and backward reactions can form a
%                   loop and therefore appear connected), but it makes the
%                   problem significantly more computationally intensive to
%                   solve (two more integer constraints per reversible reaction)
%                   (opt, default false)
%   params          parameter structure as used by getMILPParams (opt,
%                   default [])
%
%   outModel        the resulting model structure
%   deletedRxns     reactions which were deleted by the algorithm
%   metProduction   array that indicates which of the
%                   metabolites in presentMets that could be
%                   produced
%                   -2: metabolite name not found in model
%                   -1: metabolite found, but it could not be produced
%                   1: metabolite could be produced
%   fValue          objective value (sum of (the negative of)
%                   reaction scores for the included reactions and
%                   prodWeight*number of produced metabolites)
%
%   This function is the actual implementation of the algorithm. See
%   getINITModel for a higher-level function for model reconstruction. See
%   PLoS Comput Biol. 2012;8(5):e1002518 for details regarding the
%   implementation.
%
%   Usage: [outModel deletedRxns metProduction fValue]=runINIT(model,...
%           rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,...
%           noRevLoops,params)

if nargin<2
    rxnScores=zeros(numel(model.rxns),1);
end
if isempty(rxnScores)
    rxnScores=zeros(numel(model.rxns),1);
end
if nargin<3
    presentMets={};
end
if isempty(presentMets)
    presentMets={};
end
presentMets=presentMets(:);
if nargin<4
    essentialRxns={};
end
if isempty(essentialRxns)
    essentialRxns={};
end
essentialRxns=essentialRxns(:);
if nargin<5
    prodWeight=0.5;
end
if isempty(prodWeight)
    prodWeight=0.5;
end
if nargin<6
    allowExcretion=false;
end
if nargin<7
    noRevLoops=false;
end
if nargin<8
    params=[];
end

if numel(presentMets)~=numel(unique(presentMets))
    EM='Duplicate metabolite names in presentMets';
    dispEM(EM);
end

%Default is that the metabolites cannot be produced
if ~isempty(presentMets)
    metProduction=ones(numel(presentMets),1)*-2;
    presentMets=upper(presentMets);
    pmIndexes=find(ismember(presentMets,upper(model.metNames)));
    metProduction(pmIndexes)=-1; %Then set that they are at least found
else
    metProduction=[];
    pmIndexes=[];
end

%The model should be in the reversible format and all relevant exchange
%reactions should be open
if isfield(model,'unconstrained')
    EM='Exchange metabolites are still present in the model. Use simplifyModel if this is not intended';
    dispEM(EM,false);
end

%The irreversible reactions that are essential must have a flux and are
%therefore not optimized for using MILP, which reduces the problem size.
%However, reversible reactions must have a flux in one direction, so they
%have to stay in the problem. The essentiality constraint on reversible
%reactions is implemented in the same manner as for reversible reactions
%when noRevLoops==true, but with the additional constraint that C ub=-1.
%This forces one of the directions to be active.
revRxns=find(model.rev~=0);
essentialReversible=find(ismember(model.rxns(revRxns),essentialRxns));
essentialRxns=intersect(essentialRxns,model.rxns(model.rev==0));

%Convert the model to irreversible
irrevModel=convertToIrrev(model);
rxnScores=[rxnScores;rxnScores(model.rev==1)];
%These are used if noRevLoops is true
if noRevLoops==true
    forwardIndexes=find(model.rev~=0);
    backwardIndexes=(numel(model.rxns)+1:numel(irrevModel.rxns))';
else
    %Then they should only be used for essential reversible reactions
    forwardIndexes=revRxns(essentialReversible);
    backwardIndexes=essentialReversible+numel(model.rxns);
end

%Get the indexes of the essential reactions and remove them from the
%scoring vector
essentialIndex=find(ismember(irrevModel.rxns,essentialRxns));
rxnScores(essentialIndex)=[];

%Go through each of the presentMets (if they exist) and modify the S matrix
%so that each reaction which produces any of them also produces a
%corresponding fake metabolite and the opposite in the reverse direction.

%This is to deal with the fact that there is no compartment info regarding
%the presentMets. This modifies the irrevModel structure, but that is fine
%since it's the model structure that is returned.
if any(pmIndexes)
    irrevModel.metNames=upper(irrevModel.metNames);
    metsToAdd.mets=strcat({'FAKEFORPM'},num2str(pmIndexes));
    metsToAdd.metNames=metsToAdd.mets;
    metsToAdd.compartments=irrevModel.comps{1};
    
    %There is no constraints on the metabolites yet, since maybe not all of
    %them could be produced
    irrevModel=addMets(irrevModel,metsToAdd);
end

%Modify the matrix
for i=1:numel(pmIndexes)
    %Get the matching mets
    I=ismember(irrevModel.metNames,presentMets(pmIndexes(i)));
    
    %Find the reactions where any of them are used.
    [~, K, L]=find(irrevModel.S(I,:));
    
    %This ugly loop is to avoid problems if a metabolite occurs several
    %times in one reaction
    KK=unique(K);
    LL=zeros(numel(KK),1);
    for j=1:numel(KK)
        LL(j)=sum(L(K==KK(j)));
    end
    irrevModel.S(numel(irrevModel.mets)-numel(pmIndexes)+i,KK)=LL;
end

%Some nice to have numbers
nMets=numel(irrevModel.mets);
nRxns=numel(irrevModel.rxns);
nEssential=numel(essentialIndex);
nNonEssential=nRxns-nEssential;
nonEssentialIndex=setdiff(1:nRxns,essentialIndex);
S=irrevModel.S;

%Add so that each non-essential reaction produces one unit of a fake
%metabolite
temp=sparse(1:nRxns,1:nRxns,1);
temp(essentialIndex,:)=[];
S=[S;temp];

%Add another set of reactions (will be binary) which also produce these
%fake metabolites, but with a stoichiometry of 1000
temp=sparse(1:nNonEssential,1:nNonEssential,1000);
temp=[sparse(nMets,nNonEssential);temp];
S=[S temp];

%Add reactions for net-production of (real) metabolites
if prodWeight~=0
    temp=[speye(nMets-numel(pmIndexes))*-1;sparse(nNonEssential+numel(pmIndexes),nMets-numel(pmIndexes))];
    S=[S temp];
    %To keep the number of reactions added like this
    nNetProd=nMets-numel(pmIndexes);
else
    nNetProd=0;
end

%Add constraints so that reversible reactions can only be used in one
%direction. This is done by adding the fake metabolites A, B, C for each
%reversible reaction in the following manner
% forward: A + .. => ... backwards: B + ... => ... int1: C => 1000 A int2:
% C => 1000 B A ub=999.9 B ub=999.9 C lb=-1 int1 and int2 are binary
if any(forwardIndexes)
    nRevBounds=numel(forwardIndexes);
    
    %Add the A metabolites for the forward reactions and the B metabolites
    %for the reverse reactions
    I=speye(numel(irrevModel.rxns))*-1;
    temp=[I(forwardIndexes,:);I(backwardIndexes,:)];
    
    %Padding
    temp=[temp sparse(size(temp,1),size(S,2)-numel(irrevModel.rxns))];
    
    %Add the int1 & int2 reactions that produce A and B
    temp=[temp speye(nRevBounds*2)*1000];
    
    %And add that they also consume C
    temp=[temp;[sparse(nRevBounds,size(S,2)) speye(nRevBounds)*-1 speye(nRevBounds)*-1]];
    
    %Add the new reactions and metabolites
    S=[S sparse(size(S,1),nRevBounds*2)];
    S=[S;temp];
else
    nRevBounds=0;
end

%Add so that the essential reactions must have a small flux and that the
%binary ones (and net-production reactions) may have zero flux. The integer
%reactions for reversible reactions have [0 1]
prob.blx=[irrevModel.lb;zeros(nNonEssential+nNetProd+nRevBounds*2,1)];
prob.blx(essentialIndex)=max(0.1,prob.blx(essentialIndex));

%Add so that the binary ones and net-production reactions can have at the
%most flux 1.0
prob.bux=[irrevModel.ub;ones(nNonEssential+nNetProd+nRevBounds*2,1)];

%Add that the fake metabolites must be produced in a small amount and that
%the A and B metabolites for reversible reactions can be [0 999.9] and C
%metabolites [-1 0]
prob.blc=[irrevModel.b(:,1);ones(nNonEssential,1);zeros(nRevBounds*2,1);ones(nRevBounds,1)*-1];

%Add that normal metabolites can be freely excreted if
%allowExcretion==true, and that the fake ones can be excreted 1000 units at
%most. C metabolites for essential reversible reactions should have an
%upper bound of -1. If noRevLoops is false, then add this constraint for
%all the reactions instead.
if noRevLoops==true
    revUB=zeros(nRevBounds,1);
    revUB(essentialReversible)=-1;
else
    revUB=ones(nRevBounds,1)*-1;
end
if allowExcretion==true
    metUB=inf(nMets,1);
else
    metUB=irrevModel.b(:,min(size(irrevModel.b,2),2));
end
prob.buc=[metUB;ones(nNonEssential,1)*1000;ones(nRevBounds*2,1)*999.9;revUB];

%Add objective coefficients for the binary reactions. The negative is used
%since we're minimizing. The negative is taken for the prodWeight as well,
%in order to be consistent with the syntax that positive scores are good
prob.c=[zeros(nRxns,1);rxnScores;ones(nNetProd,1)*prodWeight*-1;zeros(nRevBounds*2,1)];
prob.a=S;

% adapt problem structure for cobra-style solver
prob.c=[prob.c;zeros(size(prob.a,1),1)];
prob.A=[prob.a -speye(size(prob.a,1))];
prob.b=zeros(size(prob.a,1), 1);
prob.ub=[prob.bux; prob.buc];
prob.osense=1;
prob.csense=char(zeros(size(prob.a,1),1));
prob.csense(:)='E';

%We still don't know which of the presentMets that can be produced. Go
%through them, force production, and see if the problem can be solved
for i=1:numel(pmIndexes)
    prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=1;
    prob.lb=[prob.blx; prob.blc];
    res=optimizeProb(prob,params);
    isFeasible=checkSolution(res);
    if ~isFeasible
        %Reset the constraint again
        prob.blc(numel(irrevModel.mets)-numel(pmIndexes)+i)=0;
    else
        %Metabolite produced
        metProduction(pmIndexes(i))=1;
    end
end
prob.lb=[prob.blx; prob.blc];

%Add that the binary reactions may only take integer values.
prob.vartype = repmat('C', 1, size(prob.A, 2));
allInt=[(nRxns+1):(nRxns+nNonEssential) size(S,2)-nRevBounds*2+1:size(S,2)];
prob.vartype(allInt) = 'B';

% solve problem
res=optimizeProb(prob,params);

%Problem should not be infeasible, but it is possible that the time limit
%was reached before finding any solutions.
if ~checkSolution(res)
    if strcmp(res.origStat, 'TIME_LIMIT')
        EM='Time limit reached without finding a solution. Try increasing the TimeLimit parameter.';
    else
        EM='The problem is infeasible';
    end
    dispEM(EM);
end

fValue=res.obj;

%Get all reactions used in the irreversible model
usedRxns=(nonEssentialIndex(res.full(nRxns+1:nRxns+nNonEssential)<0.1))';

%Map to reversible model IDs
usedRxns=[usedRxns(usedRxns<=numel(model.rxns));revRxns(usedRxns(usedRxns>numel(model.rxns))-numel(model.rxns))];

%Then get the ones that are not used in either direction or is essential
I=true(numel(model.rxns),1);
I(usedRxns)=false;
I(essentialIndex)=false;
deletedRxns=model.rxns(I);
outModel=removeReactions(model,I,true,true);
end
