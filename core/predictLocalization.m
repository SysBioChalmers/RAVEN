function [outModel, geneLocalization, transportStruct, scores, removedRxns]=predictLocalization(model,GSS,defaultCompartment,transportCost,maxTime,plotResults)
% predictLocalization
%   Tries to assign reactions to compartments in a manner that is in
%   agreement with localization predictors while at the same time
%   maintaining connectivity.
%
%   model                 a model structure. If the model contains several
%                         compartments they will be merged
%   GSS                   gene scoring structure as from parseScores
%   defaultCompartment    transport reactions are expressed as diffusion
%                         between the defaultCompartment and the others.
%                         This is usually the cytosol. The default
%                         compartment must have a match in
%                         GSS
%   transportCost         the cost for including a transport reaction. If this
%                         a scalar then the same cost is used for all metabolites.
%                         It can also be a vector of costs with the same dimension
%                         as model.mets. Note that negative costs will result in that
%                         transport of the metabolite is encouraged (opt, default 0.5)
%   maxTime               maximum optimization time in minutes (opt,
%                         default 15)
%   plotResults           true if the result should be plotted during the
%                         optimization (opt false)
%
%   outModel              the resulting model structure
%   geneLocalization      structure with the genes and their resulting
%                         localization
%   transportStruct       structure with the transport reactions that had
%                         to be inferred and between which compartments
%   scores                structure that contains the total score history
%                         together with the score based on gene localization
%                         and the score based on included transport reactions
%   removedRxns           cell array with the reaction ids that had to be
%                         removed in order to have a connected input model
%
%   This function requires that the starting network is connected when it's in
%   one compartment. Reactions that are unconnected are removed and saved
%   in removedRxns. Try running fillGaps to have a more connected input
%   model if there are many such reactions.
%
%   In the final model all metabolites are produced in at least one reaction.
%   This doesn't guarantee a fully functional model since there can be internal
%   loops. Transport reactions are only included as passive diffusion (A <=> B).
%
%   The score of a model is the sum of scores for all genes in their
%   assigned compartment minus the cost of all transport reactions that
%   had to be included. A gene can only be assigned to one compartment.
%   This is a simplification to keep the problem size down. The problem is
%   solved using simulated annealing.
%
%   Usage: [outModel, geneLocalization, transportStruct, score, removedRxns]=...
%       predictLocalization(model,GSS,defaultCompartment,transportCost,maxTime)
%
%   Simonas Marcisauskas, 2017-08-25
%

if nargin<4
    transportCost=ones(numel(model.mets),1)*0.5;
end
if numel(transportCost)==1
    transportCost=ones(numel(model.mets),1)*transportCost;
end
transportCost=transportCost(:);

if numel(transportCost)~=numel(model.mets)
    EM='The vector of transport costs must have the same dimension as model.mets';
    dispEM(EM,true);
end
if nargin<5
    maxTime=15;
end
if nargin<6
    plotResults=false;
end

if isfield(model,'rxnComps')
   model=rmfield(model,'rxnComps');
   EM='The model structure contains information about reaction compartmentalization. This is not supported by this function. The rxnComps field has been deleted';
   dispEM(EM,false);
end
if isfield(model,'geneComps')
   model=rmfield(model,'geneComps');
   EM='The model structure contains information about gene compartmentalization. This is not supported by this function. The geneComps field has been deleted';
   dispEM(EM,false);
end

I=ismember(defaultCompartment,GSS.compartments);
if I==false
    EM='defaultCompartment not found in GSS';
    dispEM(EM);
end

if numel(model.comps)>1
    EM='The model has several compartments. All compartments will be merged';
    dispEM(EM,false);
    model=mergeCompartments(model,true,true);
end

%***Begin formating the data structures

%Expand the model so that iso-enzymes have different reactions
model=expandModel(model);

%Identify reactions that have to be deleted because the involved mets are
%never produced. This is done in an iterative manner
removedRxns={};
%This is to keep track of which metabolites are removed in this step. It is
%needed to adjust the transport costs
originalModelMets=model.mets;
while 1
    irrevModel=convertToIrrev(model);

    I=sum(irrevModel.S>0,2);

    %Pretend that the unconstrained metabolites are made enough
    if isfield(irrevModel,'unconstrained')
        I(irrevModel.unconstrained~=0)=2;
    end
    metsToDelete=false(numel(model.mets),1);

    %This is not very neat but I loop through each metabolite and check whether
    %it can be produced (without using only one isolated reversible reaction)
    for i=1:numel(irrevModel.mets)
        %If something can be made in two reactions then everything is fine. If i
        %can be made in one reaction it's fine unless it's through an isolated
        %reversible reaction (which can act as a mini loop)
        if I(i)<2
            if I(i)==1
                %Find the reaction where this metabolite is produced
                [~, J]=find(irrevModel.S(i,:)>0);

                %Check the metabolites that are consumed in this reaction. The
                %problem is if any of them is only produced in the opposite
                %reversible reaction
                K=irrevModel.S(:,J)<0;
                check=find(K & I<=1);

                for j=1:numel(check)
                    %Find the reactions where it participates
                    [~, L]=find(irrevModel.S(check(j),:)>0);

                    if ~isempty(L)
                        rxn=irrevModel.rxns(J);
                        rxnRev=irrevModel.rxns(L);
                        if strcmp(strrep(rxn,'_REV',''),strrep(rxnRev,'_REV',''))
                           metsToDelete(i)=true;
                        end
                    else
                        %If the metabolite was never produced then do nothing and deal with
                        %it when the loop gets there :)
                        continue;
                    end
                end
            else
                %Not made anywhere
                metsToDelete(i)=true;
            end
        end
    end

    if any(metsToDelete)
        %Delete any reactions involving any of the metsToDelete
        [~, I]=find(model.S(metsToDelete,:));
        removedRxns=[removedRxns;model.rxns(I)];
        model=removeReactions(model,I,true,true);
    else
        %All bad reactions deleted
        break;
    end
end

%Adjust the transport costs
transportCost=transportCost(ismember(originalModelMets,model.mets));

%Assign fake genes to reactions without genes. This is just to make things
%easier later on
I=find(sum(model.rxnGeneMat,2)==0);
for i=1:numel(I)
    model.genes=[model.genes;['&&FAKE&&' num2str(i)]];
    if isfield(model,'geneShortNames')
        model.geneShortNames=[model.geneShortNames;{''}];
    end
    if isfield(model,'geneMiriams')
        model.geneMiriams=[model.geneMiriams;{[]}];
    end
    if isfield(model,'geneFrom')
        model.geneFrom=[model.geneFrom;{{'FAKE'}}];
    end
    model.rxnGeneMat(I(i),numel(model.genes))=1;
    model.grRules{I(i)}='';
end

%Update the GSS. All genes, fake or real, for which
%there is no evidence gets a score 0.5 in all compartments. Also just to
%make it easier further on
I=setdiff(model.genes,GSS.genes);
GSS.genes=[GSS.genes;I];
GSS.scores=[GSS.scores;ones(numel(I),numel(GSS.compartments))*0.5];

%Gene complexes should be moved together in order to be biologically
%relevant. The average score for the genes is used for each compartment.
%This is done by changing the model so that gene complexes are used as a
%single gene name and then a score is calculated for that "gene".

%Only "and"-relationships exist after expandModel
genes=unique(model.grRules);
nGenes=strrep(genes,'(','');
nGenes=strrep(nGenes,')','');
%nGenes=strrep(nGenes,' and ','_and_');
complexes=setdiff(nGenes,model.genes);
if ~isempty(complexes)
    if isempty(complexes{1}) %Empty grRules also come up here
        complexes(1)=[];
    end
end
cScores=zeros(numel(complexes),numel(GSS.compartments));
for i=1:numel(complexes)
    genesInComplex=regexp(complexes{i},' and ','split');

    %Find these genes in GSS
    [I, J]=ismember(genesInComplex,GSS.genes);

    if any(I)
        %Get the average of the genes that were found.
        mScores=mean(GSS.scores(J(I),:));

        %And add 0.5 for the genes that were not found in order to be
        %consistent with non-complexes
        mScores=(mScores.*sum(I)+(numel(genesInComplex)-sum(I))*0.5)/numel(genesInComplex);
    else
        EM=['Could not parse grRule "' complexes{i} '". Assigning score 0.0 in all compartments'];
        dispEM(EM,false);
        mScores=ones(1,numel(genesInComplex))*0.5;
    end
    cScores(i,:)=mScores;

    %Add this complex as a new gene
    model.genes=[model.genes;complexes{i}];
    if isfield(model,'geneMiriams')
        model.geneMiriams=[model.geneMiriams;{[]}];
    end
    if isfield(model,'geneShortNames')
        model.geneShortNames=[model.geneShortNames;{''}];
    end
    if isfield(model,'geneFrom')
        model.geneFrom=[model.geneFrom;{'COMPLEX'}];
    end
    %Find the reactions which had the original complex and change them to
    %use the new "gene"
    I=ismember(model.grRules,['(' complexes{i} ')']);

    %Should check more carefully if there can be an error here
    if ~isempty(I)
        model.rxnGeneMat(I,:)=0; %Ok since we have split on "or"
        model.rxnGeneMat(I,numel(model.genes))=1;
    end
end

%Add the new "genes"
GSS.genes=[GSS.genes;complexes];
GSS.scores=[GSS.scores;cScores];

%After merging the complexes it could happen that there are genes that are
%no longer in use. Delete such genes
model=removeReactions(model,{},false,true);

%Exchange reactions, defined as involving an unconstrained metabolite, are
%special in that they have to stay in the defaultCompartment. This means
%that uptake/excretion of metabolites is always via the default
%compartment. This is a small simplification, but should be valid in most
%cases
[~, I]=getExchangeRxns(model);

%It will be easier later on if the same place. Put them in the beginning
J=1:numel(model.rxns);
J(I)=[];
model=permuteModel(model,[I;J'],'rxns');

%Number of exchange reactions
nER=numel(I);

%Also put the exchange metabolites in the beginning
if isfield(model,'unconstrained')
    I=find(model.unconstrained);
    J=1:numel(model.mets);
    J(I)=[];
    model=permuteModel(model,[I;J'],'mets');
    %Also reorder the transport costs
    transportCost=transportCost([I;J']);
    %Number of exchange metabolites
    nEM=numel(I);
else
    nEM=0;
end

%There is no point of having genes for exchange reactions, so delete them.
%Also to make computations easier.
model.rxnGeneMat(1:nER,:)=0;
model.grRules(1:nER)={''};

%Remove unused genes
model=removeReactions(model,{},false,true);

%Remove genes with no match to the model and reorder so that the genes are
%in the same order as model.genes. Since we have already added fake genes
%so that all genes in model exist in GSS it's fine to do
%like this.
[~, J]=ismember(model.genes,GSS.genes);
GSS.genes=model.genes;
GSS.scores=GSS.scores(J,:);

%Reorder the GSS so that the first index corresponds to the
%default compartment
[~, J]=ismember(defaultCompartment,GSS.compartments);
reorder=1:numel(GSS.compartments);
reorder(J)=[];
reorder=[J reorder];
GSS.scores=GSS.scores(:,reorder);
GSS.compartments=GSS.compartments(reorder);

%Since we're only looking at whether the metabolites can be synthesized, we
%don't have to care about the stoichiometry. Change to -1/1 to simplify
%later. Keep the S matrix for later though.
oldS=model.S;
model.S(model.S>0)=1;
model.S(model.S<0)=-1;

%Here I do a bit of a trick. Since I don't want to calculate which
%reactions are reversible all the time, I let reversible reactions have the
%coefficients -10/10 instead of -1/1
model.S(:,model.rev==1)=model.S(:,model.rev==1).*10;

%***Begin problem formulation

%Some numbers that are good to have
nRxns=numel(model.rxns)-nER; %Excluding exchange rxns
nMets=numel(model.mets)-nEM; %Excluding exchange mets
nGenes=numel(model.genes);
nComps=numel(GSS.compartments);

%Create a big stoichiometric matrix that will be the current model. In
%order to have faster simulations the maximal model size is declared and
%reactions are then moved within it.

%First the original model (with the first nE being exchange rxns), then
%reserve space for number of rxns minus exchange rxns for each non-default
%compartment, then transport reactions for all non-exchange mets between
%the default compartment and all others.
%NOTE: Kept eye()*0 since eye() can be used to include all transport from
%the beginning
s=repmat(eye(nMets)*0,1,nComps-1);
s=[zeros(numel(model.mets)-nMets,size(s,2));s];
S=[model.S sparse(numel(model.mets),nRxns*(nComps-1)) s];
s=[sparse(nMets*(nComps-1),numel(model.rxns)+nRxns*(nComps-1)) eye(nMets*(nComps-1))*0];
S=[S;s];

%Also replicate the transport costs
transportCost=[transportCost(1:nEM);repmat(transportCost(nEM+1:end),nComps,1)];

%Create a binary matrix that says where the genes are in the current
%solution
g2c=false(nGenes,nComps);
%All genes start in the default compartment
g2c(:,1)=true;

%Start of main optimization loop
tic;
bestScore=-inf;
bestS=[];
bestg2c=[];

%Temp for testing
plotScore=[];
nTrans=[];
totScore=[];
minScore=sum(min(GSS.scores,[],2));
maxScore=sum(max(GSS.scores,[],2));

while toc<maxTime*60
   %Pick a random gene, weighted by it's current score minus the best score
   %for that gene (often 1.0, but can be 0.5 for no genes or average for complexes.
   %Genes with bad fits are more likely to be moved. This formulation never
   %moves a gene from its best compartment. Therefore a small uniform
   %weight is added.
   [I, J]=find(g2c);
   geneToMove=randsample(nGenes,1,true,max(GSS.scores(I,:),[],2)-GSS.scores(sub2ind(size(g2c),I,J))+0.1);

   %Sample among possible compartments to move to. Add a larger weight to
   %even out the odds a little. Also a way of getting rid of loops where
   %the same set of genes are moved back and forth several times.
   toComp=randsample(nComps,1,true,GSS.scores(geneToMove,:)+0.2);

   %Check that it moves to a new compartment
   if toComp==find(g2c(geneToMove,:))
       continue;
   end

   %Moves the gene
   [newS, newg2c]=moveGene(S,model,g2c,geneToMove,toComp,nRxns,nMets);

   %Tries to connect the network. If this was not possible in 10
   %iterations, then abort. If more than 20 modifications were needed then
   %it's unlikely that it will be a lower score
   wasConnected=false;
   for j=1:10
       %Find the metabolites that are now unconnected
       unconnected=findUnconnected(newS,nEM);

       %Continue if there are still unconnected
       if any(unconnected)
           %For each gene find out how many of these could be connected if
           %the gene was moved and how many would be disconnected by moving
           %that gene
           [geneIndex, moveTo, deltaConnected, deltaScore]=selectGenes(newS,nEM,nMets,nER,nRxns,model,unconnected,g2c,GSS);

           %Score which gene would be the best to move. The highest
           %deltaScore is 1.0. I want it to be possible to move a gene from
           %worst to best compartment even if it disconnects, say, 1.5 more
           %metabolites.
           [score, I]=max(1.5*deltaScore+deltaConnected);

           %Checks if it has to add a transport or if there is a gene that
           %could be moved order to have a more connected network
           hasToAddTransport=true;
           if ~isempty(deltaConnected)
              if score>0
                  hasToAddTransport=false;
              end
           end

           %If it is possible to move any gene in order to have a more
           %connected network, then move the best one
           if hasToAddTransport==false;
                [newS, newg2c]=moveGene(newS,model,g2c,geneIndex(I),moveTo(I),nRxns,nMets);
           else
                %Choose a random unconnected metabolite that should be
                %connected
                transMet=unconnected(randsample(numel(unconnected),1));

                %First get where the metabolite is now
                comps=ceil((transMet-nEM)/((size(S,1)-nEM)/nComps));

                %Find the corresponding metabolite index if it were in the
                %default compartment
                dcIndex=transMet-(comps-1)*nMets;

                %Then get the indexes of that metabolite in all
                %compartments
                allIndexes=dcIndex;
                for k=1:nComps-1
                   allIndexes=[allIndexes;dcIndex+nMets*k];
                end

                %It could be that some of these aren't used in any
                %reaction. Get only the ones which are
                I=sum(newS(allIndexes,:)~=0,2)>0;

                %Then get the ones that are used but not in unconnected.
                %These are metabolites that could potentially be
                %transported to connect transMet
                connectedUsed=setdiff(allIndexes(I),unconnected);

                %I think this is an error but I leave it for now. It seems
                %to happen if nothing can be connected in one step
                if isempty(connectedUsed)
                   break;
                end

                %If transMet is in the default compartment then everything
                %is fine, just connect it to a random one
                if transMet==dcIndex
                	newS=addTransport(newS,nRxns,nER,nMets,nEM,nComps,transMet,connectedUsed(randsample(numel(connectedUsed),1)));
                else
                    %If one of the connectedUsed is in the default
                    %compartment then connect to that one
                    I=connectedUsed(connectedUsed<(nMets+nEM));
                    if any(I)
                    	newS=addTransport(newS,nRxns,nER,nMets,nEM,nComps,transMet,I(randsample(numel(I),1)));
                    else
                        %This is if the only way to connect it is by adding
                        %two transport reactions, going via the default
                        %compartment
                        break;
                    end
                end
           end
       else
           wasConnected=true;
           break;
       end
   end

   %If the network was connected in a new way, it is possible that some
   %transport reactions are no longer needed. They should be removed
   if wasConnected==true
        %These are the metabolites that are being transported
        activeTransport=find(sum(newS(:,nER+nRxns*nComps+1:end),2));

        %Get the metabolites that are unconnected if transport wasn't used
        unconnected=findUnconnected(newS(:,1:nER+nRxns*nComps),nEM);

        %Find the transport reactions that are not needed and delete them
        I=setdiff(activeTransport,unconnected);

        %Since both metabolites in a transport rxns must be connected for
        %the reaction to be deleted, the sum over the colums should be 4.
        newS(:,find(sum(newS(I,nER+nRxns*nComps+1:end))==4)+nER+nRxns*nComps)=0;

        %Score the solution and determine whether to keep it as a new solution
        [score, geneScore, trCost]=scoreModel(newS,newg2c,GSS,transportCost);

        %If it was the best solution so far, keep it
        if score>bestScore
            bestScore=score;
            bestS=newS;
            bestg2c=newg2c;
        end

        %This should not be steepest descent later
        if score>=bestScore% || exp((score-bestScore)*7)>rand()
            plotScore=[plotScore;geneScore];
            nTrans=[nTrans;trCost];
            totScore=[totScore;score];
            S=newS;
            g2c=newg2c;

            if plotResults==true
                subplot(3,2,1);
                spy(S);
                subplot(3,2,2);
                plot(plotScore,'r');
                xlabel('Gene score');
                subplot(3,2,3);
                plot((plotScore-minScore)/(maxScore-minScore),'r');
                xlabel('Gene score relative to predictions');
                subplot(3,2,4);
                plot(nTrans,'g');
                xlabel('Transport cost');
                subplot(3,2,5);
                plot(totScore,'b');
                xlabel('Total score');
                subplot(3,2,6);
                pause(0.2);
            end
        end
   end
end
scores.totScore=score;
scores.geneScore=geneScore;
scores.transCost=trCost;

%Find which metabolites are transported and to where
[I, J]=find(bestS(nEM+1:nEM+nMets,end-nMets*(nComps-1)+1:end));
J=ceil(J/nMets+1);
transportStruct.mets=model.metNames(I+nEM);
transportStruct.toComp=GSS.compartments(J);

[I, J]=find(bestg2c);
geneLocalization.genes=GSS.genes(I);
geneLocalization.comps=GSS.compartments(J);

%Resort the gene names
[~, I]=sort(geneLocalization.genes);
geneLocalization.genes=geneLocalization.genes(I);
geneLocalization.comps=geneLocalization.comps(I);

%Remove the fake genes
I=strncmp('&&FAKE&&',geneLocalization.genes,8);
geneLocalization.genes(I)=[];
geneLocalization.comps(I)=[];

%Put together the model. This is done by first duplicating the S matrix into the
%different compartments. Then the transport reactions are added based on
%transportStruct. By now model.S should have the same size as the S matrix
%used in the optimization, but with conserved stoichiometry. In the final
%step all reactions and metabolites that aren't used in the S matrix from the optimization
%are deleted from the model.
outModel=model;
outModel.S=oldS;

%This is the S matrix without exchange rxns or metabolites
copyPart=outModel.S(nEM+1:end,nER+1:end);

%Replicate to give the rxnGeneMat for the full system
copyRxnGeneMat=outModel.rxnGeneMat(nER+1:end,:);
outModel.rxnGeneMat=[outModel.rxnGeneMat;repmat(copyRxnGeneMat,nComps-1,1)];

%First fix the compartments. The model is already ordered with the exchange
%metabolites first. The original model may contain one or two compartments,
%depending on whether any exchange metabolites are defined.
nStartComps=numel(outModel.comps);
if nStartComps==1
   outModel.comps={'1'};
   outModel.compNames=GSS.compartments(1);
else
    if model.metComps(1)==1
        outModel.compNames(1)=GSS.compartments(1);
    else
        outModel.compNames(2)=GSS.compartments(1);
    end
end
outModel.compNames=[outModel.compNames;GSS.compartments(2:end)];

%Ugly little loop
for i=1:numel(GSS.compartments)-1
    outModel.comps=[outModel.comps;num2str(numel(outModel.comps)+1)];
end
%This information is not known from the data, so empty fields are added
outModel.compOutside=cell(numel(outModel.comps),1);
outModel.compOutside(:)={''};

for i=1:nComps-1
    outModel.S=[outModel.S sparse(size(outModel.S,1),nRxns)];
    outModel.S=[outModel.S;[sparse(nMets,nRxns*i+nER) copyPart]];
    outModel.rxns=[outModel.rxns;strcat(outModel.rxns(nER+1:nER+nRxns),'_',GSS.compartments{i+1})];
    outModel.rxnNames=[outModel.rxnNames;strcat(outModel.rxnNames(nER+1:nER+nRxns),' (',GSS.compartments{i+1},')')];
    outModel.lb=[outModel.lb;outModel.lb(nER+1:nER+nRxns)];
    outModel.ub=[outModel.ub;outModel.ub(nER+1:nER+nRxns)];
    outModel.rev=[outModel.rev;outModel.rev(nER+1:nER+nRxns)];
    outModel.c=[outModel.c;outModel.c(nER+1:nER+nRxns)];
    if isfield(outModel,'grRules')
        outModel.grRules=[outModel.grRules;outModel.grRules(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'subSystems')
        outModel.subSystems=[outModel.subSystems;outModel.subSystems(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'eccodes')
        outModel.eccodes=[outModel.eccodes;outModel.eccodes(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'rxnFrom')
        outModel.rxnFrom=[outModel.rxnFrom;outModel.rxnFrom(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'rxnMiriams')
        outModel.rxnMiriams=[outModel.rxnMiriams;outModel.rxnMiriams(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'rxnNotes')
        outModel.rxnNotes=[outModel.rxnNotes;outModel.rxnNotes(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'rxnReferences')
        outModel.rxnReferences=[outModel.rxnReferences;outModel.rxnReferences(nER+1:nER+nRxns)];
    end
    if isfield(outModel,'rxnConfidenceScores')
        outModel.rxnConfidenceScores=[outModel.rxnConfidenceScores;outModel.rxnConfidenceScores(nER+1:nER+nRxns)];
    end
    outModel.mets=[outModel.mets;strcat(outModel.mets(nEM+1:nEM+nMets),'_',GSS.compartments{i+1})];
    outModel.metNames=[outModel.metNames;outModel.metNames(nEM+1:nEM+nMets)];
    outModel.b=[outModel.b;outModel.b(nEM+1:nEM+nMets,:)];
    I=ones(nMets,1)*nStartComps+i;
    outModel.metComps=[outModel.metComps;I];
    if isfield(outModel,'inchis')
        outModel.inchis=[outModel.inchis;outModel.inchis(nEM+1:nEM+nMets)];
    end
    if isfield(outModel,'unconstrained')
        outModel.unconstrained=[outModel.unconstrained;outModel.unconstrained(nEM+1:nEM+nMets)];
    end
    if isfield(outModel,'metMiriams')
        outModel.metMiriams=[outModel.metMiriams;outModel.metMiriams(nEM+1:nEM+nMets)];
    end
    if isfield(outModel,'metFormulas')
        outModel.metFormulas=[outModel.metFormulas;outModel.metFormulas(nEM+1:nEM+nMets)];
    end
    if isfield(outModel,'metFrom')
        outModel.metFrom=[outModel.metFrom;outModel.metFrom(nEM+1:nEM+nMets)];
    end
    if isfield(outModel,'metCharges')
        outModel.metCharges=[outModel.metCharges;outModel.metCharges(nEM+1:nEM+nMets)];
    end
end

%Add the transport reactions
transS=bestS(:,numel(outModel.rxns)+1:end);
J=sum(transS)>0; %Active rxns

%Transport reactions are written in a different way compared to a "real"
%stoichimetric matrix. This is to fix that
transS(transS~=0)=1;
transS(1:nEM+nMets,:)=transS(1:nEM+nMets,:)*-1;
I=find(sum(transS>0,2));
nTransRxns=numel(I);
outModel.S=[outModel.S transS(:,J)];
filler=ones(nTransRxns,1);
outModel.lb=[outModel.lb;filler*-1000];
outModel.ub=[outModel.ub;filler*1000];
outModel.rev=[outModel.rev;filler];
outModel.c=[outModel.c;filler*0];
outModel.rxnGeneMat=[outModel.rxnGeneMat;sparse(nTransRxns,numel(outModel.genes))];

for i=1:numel(I)
    outModel.rxns=[outModel.rxns;strcat('transport',num2str(i))];
    outModel.rxnNames=[outModel.rxnNames;['Transport of ',outModel.metNames{I(i)}]];
    if isfield(outModel,'grRules')
        outModel.grRules=[outModel.grRules;{''}];
    end
    if isfield(outModel,'rxnMiriams')
        outModel.rxnMiriams=[outModel.rxnMiriams;{[]}];
    end
    if isfield(outModel,'subSystems')
        outModel.subSystems=[outModel.subSystems;'Inferred transport reactions'];
    end
    if isfield(outModel,'eccodes')
        outModel.eccodes=[outModel.eccodes;{''}];
    end
    if isfield(outModel,'rxnFrom')
        outModel.rxnFrom=[outModel.rxnFrom;{''}];
    end
    if isfield(outModel,'rxnNotes')
        outModel.rxnNotes=[outModel.rxnNotes;{''}];
    end
    if isfield(outModel,'rxnReferences')
        outModel.rxnReferences=[outModel.rxnReferences;{''}];
    end
    if isfield(outModel,'rxnConfidenceScores')
        outModel.rxnConfidenceScores=[outModel.rxnConfidenceScores;{''}];
    end
end

%Then remove all reactions and metabolites that aren't used in the final
%solution from the optimization
[~, J]=find(bestS(:,1:nER+nComps*nRxns));
K=true(numel(outModel.rxns),1);
K(J)=false;
K(end-nTransRxns+1:end)=false;
outModel=removeReactions(outModel,K,true);

%Remove all fake genes
I=strncmp('&&FAKE&&',outModel.genes,8);
outModel.genes(I)=[];
if isfield(outModel,'geneMiriams')
    outModel.geneMiriams(I)=[];
end
if isfield(outModel,'geneShortNames')
    outModel.geneShortNames(I)=[];
end
outModel.rxnGeneMat(:,I)=[];
end

%Moves a gene and all associated reactions from one compartment to another
function [S, g2c]=moveGene(S,model,g2c,geneToMove,toComp,nRxns,nMets)
    %Find the current compartment and update to the new one
    currentComp=find(g2c(geneToMove,:));
    g2c(geneToMove,:)=false;
    g2c(geneToMove,toComp)=true;

    %Find the reactions in the original model that the gene controls
    [I, ~]=find(model.rxnGeneMat(:,geneToMove));

    %Calculate their current positions in the S matrix
    oldRxns=I+(currentComp-1)*nRxns;

    %And their new positions
    newRxns=I+(toComp-1)*nRxns;

    %The metabolite ids also have to be changed in order to match the new
    %compartment
    metChange=nMets*(toComp-currentComp);

    %Update the reactions
    [I, J, K]=find(S(:,oldRxns));
    I=I+metChange;

    %Move the reactions
    S(:,oldRxns)=0;
    S(sub2ind(size(S),I,newRxns(J)))=K;
end

%Finds which metabolites are unconnected, in the sense that they are never
%a product or only a product in a reversible reaction where one reactant is
%only a product in the opposite direction of that reaction. This function
%ignores exchange metabolites. Returns a vector of metabolite indexes.
%metsToCheck is an array of metabolite indexes to check for connectivity.
%If not supplied then all metabolites are checked
function unconnected=findUnconnected(S,nEM,metsToCheck)
    if nargin>2
        %Do this by deleting everything from the network that is not in
        %metsToCheck and that is not exchange metabolites
        I=false(size(S,1),1);
        I(1:nEM)=true;
        I(metsToCheck)=true;
        S=S(I,:);
    end

    em=false(size(S,1),1);
    em(1:nEM)=true;

    %Construct a matrix in which the reversible reactions are inverted
    I=sum(S>2,1) | sum(S>2,1);
    revS=S;
    revS(:,I)=revS(:,I)*-1;

    %First calculate the ones that are ok
    %Produced in 2 rxns, is exchange, is not used at all, is produced in
    %non-reversible, involved in more than 1 reversible reactions
    connected=sum(S>0,2)>1 | em | sum(S~=0,2)==0 | sum(S(:,~I)>0,2)>0 | sum(S(:,I)~=0,2)>1;

    %Then get the ones that are unconnected because they are never produced
    unconnected=sum(S>0 | revS>0,2)==0 & connected==false;

    %Then get the ones that are potentially unconnected
    maybeUnconnected=~connected & ~unconnected;
    %maybeUnconnected=find(maybeUnconnectedS);

    %The metabolites in maybeUnconnected are involved in one reversible reaction and
    %not produced in any other reaction.
    %This means that the reactions which have at least one met in
    %maybeUnconnected as reactant and one as product are unconnected. The
    %metabolites in maybeUnconnected that are present in those reactions
    %are then dead ends
    deadRxns=any(S(maybeUnconnected,:)>0) & any(S(maybeUnconnected,:)<0);

    %Get the mets involved in any of those reactions
    problematic=any(S(:,deadRxns)~=0,2);

    %If any of these are in the maybeUnconnected list then the
    %metabolite is unconnected
    unconnected(problematic & maybeUnconnected)=true;

    %Map back to metsToCheck
    if nargin>2
        unconnected=metsToCheck(unconnected(nEM+1:end));
    else
        unconnected=find(unconnected);
    end
end

%Given a set of unconnected metabolites, this function tries to move each
%gene that could connect any of them, calculates the number of newly connected
%metabolites minus the number of newly disconnected metabolites. As some metabolites
%are very connected, only 25 genes are checked. Genes that have a low score
%in their current compartment are more likely to be moved.
function [geneIndex, moveTo, deltaConnected, deltaScore]=selectGenes(S,nEM,nMets,nER,nRxns,model,unconnected,g2c,GSS)
    %If moveTo is 0 then the gene can't connect any of the metabolites
    moveTo=zeros(numel(model.genes),1);
    deltaConnected=zeros(numel(model.genes),1);

    %First get where the metabolites are now
    nComps=size(g2c,2);
    comps=ceil((unconnected-nEM)/((size(S,1)-nEM)/nComps));

	%Find the corresponding metabolite indexes if they all were in the
	%default compartment
    dcIndexes=unique(unconnected-(comps-1)*nMets);

    %Then find them if they were in any other compartment
    allIndexes=dcIndexes;
    for i=1:nComps-1
       allIndexes=[allIndexes;dcIndexes+nMets*i];
    end

    %Also check which reversible reactions that could be used
    I=sum(S>2,1) | sum(S>2,1);
    revS=S;
    revS(:,I)=revS(:,I)*-1;

    %Find all reactions that could make any of the unconnected metabolites
    %in some other compartment
    newMets=setdiff(allIndexes,unconnected);
    [~, potential]=find(S(newMets,:)>0 | revS(newMets,:)>0);
    potential(potential<=nER | potential>nER+nRxns*nComps)=[]; %No exchange rxns or transport rxns

    %Map J to the real metabolic reactions in model
    rxnComps=ceil((potential-nER)/(nRxns));

	%Find the corresponding reaction indexes if they all were in the
	%default compartment
    dcRxnIndexes=potential-(rxnComps-1)*nRxns;

    %Get the genes for those reactions
    genes=find(sum(model.rxnGeneMat(dcRxnIndexes,:)>0,1));

    %For some cases there can be very many reactions to connect something.
    %This is in particular true in the beginning of the optimization if,
    %say, ATP is unconnected. Therefore limit the number of genes to be
    %checked to 25. Weigh so that genes with bad scores in their current
    %compartment are more likely to be moved

    %Get scores for these genes
    [~, J]=find(g2c(genes,:));

    %Add a small weight so that genes in their best compartment could be
    %moved as well
    geneScores=GSS.scores(sub2ind(size(g2c),genes(:),J));
    modGeneScores=1.1-geneScores;
    if numel(genes)>25
        rGenes=genes(randsample(numel(genes),min(numel(genes),25),true,modGeneScores));

        %The sampling with weights could give duplicates
        rGenes=unique(rGenes);

        %Reorder the geneScores to match
        [~, I]=ismember(rGenes,genes);
        geneScores=geneScores(I);
        genes=rGenes;
    end
    for i=1:numel(genes)
        %Since we are moving one gene at a time, only metabolites involved
        %in any of the reactions for that gene can become unconnected. We
        %get them so speed up the algorithm.
        %First get all involved reactions in the default compartment
        rxns=find(model.rxnGeneMat(:,genes(i)));

        %Then get their mets
        mets=find(sum(model.S(:,rxns)~=0,2)>0);

        %Then get their indexes in all compartments
        allIndexes=mets;
        for j=1:nComps-1
           allIndexes=[allIndexes;mets+nMets*j];
        end

        %Check which of the unconnected metabolites that these
        %reactions correspond to. This could have been done earlier,
        %but it's fast. I skip the reversibility check because it's
        %unlikely to be an issue here. Worst case is that the gene is
        %tested once to much
        [I, ~]=find(model.S(:,rxns));
        moveToComps=unique(comps(ismember(dcIndexes,I)));

        %Try to move the gene to each of the compartments
        bestMove=-inf;
        bestComp=[];
        for j=1:numel(moveToComps)
            newS=moveGene(S,model,g2c,genes(i),moveToComps(j),nRxns,nMets);

            %Check how many metabolites that are unconnected after moving
            %the gene
            dConnected=numel(unconnected)-numel(findUnconnected(newS,nEM,[allIndexes;unconnected]));
            if dConnected>bestMove
                bestMove=dConnected;
                bestComp=moveToComps(j);
            end
        end

        %Add the difference in connectivity and where the genes should
        %be moved
        moveTo(genes(i))=bestComp;
        deltaConnected(genes(i))=bestMove;
    end

    %Finish up
    geneIndex=genes(:);
    moveTo=moveTo(geneIndex);
    deltaConnected=deltaConnected(geneIndex);
    deltaScore=GSS.scores(sub2ind(size(g2c),geneIndex(:),moveTo))-geneScores;
end

%Small function to add a transport reactions between two metabolites.
%Transport reactions are written as having a coefficient 2.0 for both
%reactant and product. This is not a "real" reaction, but since all normal
%reaction have coefficient -1/1 or -10/10 it's a compact way of writing it.
function S=addTransport(S,nRxns,nER,nMets,nEM,nComps,metA,metB)
    mets=[metA;metB];
    %Find the current compartments for the metabolites
    comps=ceil((mets-nEM)/((size(S,1)-nEM)/nComps));

    if sum(comps==1)~=1
    	EM='Tried to create a transport reaction from a non-default compartment';
      dispEM(EM);
    end

    %Calculate the reaction index
    rIndex=(nER+nRxns*nComps)+mets(comps~=1)-nEM-nMets;

    S(mets,rIndex)=2;
end

%Scores a network based on the localization of the genes and the number of
%transporter reactions used.
function [score, geneScore, transportCost]=scoreModel(S,g2c,GSS,transportCost)
    [I, J]=find(g2c);
    geneScore=sum(GSS.scores(sub2ind(size(g2c),I,J)));
    [I, ~]=find(S==2);
    I=unique(I);
    transportCost=sum(transportCost(I));
    score=geneScore-transportCost;
end
