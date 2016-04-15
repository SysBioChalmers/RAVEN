%This exercise is about creating a model from KEGG, based on protein sequences
%in a FASTA file, and doing some functionality checks on the model. The example
%case is for the yeast Saccharomyces cerevisiae. This tutorial is more of a
%showcase than the other three, and it's main purpose is to serve as a
%scaffold if you would like to reconstruct a GEM for your own organism.

%Start by downloading trained Hidden Markov Models for eukaryots from the
%RAVEN toolbox homepage. Unzip the content to a folder. Read the help for
%getKEGGModelForOrganism for details regarding the folders and what they
%are for.

%This creates a model for S. cerevisiae. The parameters are set to exclude
%general or unclear reactions and reactions with undefined stoichiometry.
%Type "help getKEGGModelForOrganism" to see what the different parameters are for.
%This process takes from several hours to a couple of days, depending on
%your hardware. The function is intended to be run on a cluster, so you can
%speed up the process by running it on several computers simultaneously.
model=getKEGGModelForOrganism('sce','sce.fa','c:\input\eukaryota','c:\output',false,false,false,10^-30,0.8,0.3,-1);

%As you can see the resulting model contains (around) 1081 reactions,
%1182 metabolites and 809 genes. Small variations is possible since it is
%an heuristic algorithm.
model


%A first control is that your model shouldn't be able to produce any metabolites
%without uptake of some metabolites. This commonly happens when metabolites
%have a different meaning in different reactions. The best way to find such
%reactions is to run makeSomething and analyze the resulting solution for
%bad reactions. An automated approach is to use removeBadRxns, which tries
%to do the same thing in an automated manner. Please type "help removeBadRxns"
%for details.
[newModel removedRxns]=removeBadRxns(model);

%You will see an error about that H+ can be made even if no reactions were
%unbalanced. Protons are particularly problematic since it's rather
%arbitary at which pH the formulas are written for. For the purpose of this
%analysis we can ignore protons and try to fix it later.
[newModel removedRxns]=removeBadRxns(model,1,{'H+'},true);

%Only one reaction was removed because it enabled the model to produce
%something from nothing. Since there were so few, it might be worthwhile to look
%into this in more detail.
removedRxns

%If you look it up in KEGG you will find that it is a general polymer reaction.
%You might want to look at the flux distributions more in detail to try to
%find out if there is any better alternative to delete. Use makeSomething
%to do this
[fluxes metabolite]=makeSomething(model,{'H+'},true);
model.metNames(metabolite)

%The model could produce H2O using the following reactions
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n')

%That resulted in quite a lot of fluxes to look through. Maybe it's easier
%if we exclude the ones that were elementally balanced.
balanceStructure=getElementalBalance(model);

%We know that hydrogen balancing is a bit tricky, so let's only look at the
%ones unbalanced for oxygen (since water was produced)
goodOnes=balanceStructure.leftComp(:,6)==balanceStructure.rightComp(:,6);
printFluxes(removeRxns(model,goodOnes), fluxes(~goodOnes), false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n')

%That didn't really help, the only unbalanced reaction was the one that was
%identified before. Let's print all fluxes involving amylose or starch
%instead.
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n',{'Amylose';'Starch'});

%This shows us two general and contradicting reactions. The first one says
%that amylose and starch are interconvertible, the second one that you can
%cleave of a glucose unit from starch to form amylose. This type of general
%reactions are problematic and should be fixed manually. We therefore
%choose to trust removeBadRxns and delete R02110.
model=removeRxns(model,'R02110');

%The model can no longer make something from nothing. Can it consume
%something without any output?
[solution metabolite]=consumeSomething(model,{'H+'},true);
model.metNames(metabolite)

%Nope, so that was good. Let's add some uptakes and see what it can
%produce.
[I J]=ismember({'D-Glucose';'H2O';'Orthophosphate';'Oxygen';'NH3';'Sulfate'},model.metNames);
[model addedRxns]=addExchangeRxns(model,'in',J);

%Check which metabolites can be produced given these uptakes. The
%canProduce function allows for output of all metabolites. This won't
%happen in the real cell, but it is very useful for functionality testing
%of the model. Once it is functional you can add excretion reactions based
%on evidence.
I=canProduce(model);

%You can see that 31% of the metabolites could be synthesized. It is not
%directly clear whether this is a high or low number, many metabolites
%should not be possible to synthesize from those simple precursors.
sum(I)/numel(model.mets)

%We can try to fill gaps using the full KEGG model to see if that gives a
%significantly higher number
keggModel=getModelFromKEGG([],false,false,false);

%The KEGG model is associated to some 800000 genes. They won't be used for
%this, so we remove them to make this a little faster
keggModel=rmfield(keggModel,'genes');
keggModel=rmfield(keggModel,'rxnGeneMat');

%We've already seen that there are some unbalanced reactions in KEGG. Let's
%only use the balanced ones for the gap filling
balanceStructure=getElementalBalance(keggModel);
keggModel=removeRxns(keggModel,balanceStructure.balanceStatus~=1,true,true);

%fillGaps with these settings will try to include reactions in order to
%have flux through all reactions in the model. There are other settings as well.
%The first flag says that production of all metabolites should be allowed.
params.relGap=0.6; %Lower number for a more exhaustive search
params.printReport=true;
[newConnected cannotConnect addedRxns newModel exitFlag]=fillGaps(model,keggModel,true,false,false,[],params);

%We see that we could connect 48 reactions (newConnected) by including 54
%reactions from the KEGG model (addedRxns). Those should of course be
%checked manually to see that they exist in yeast, but let's say that we
%have done so now.

%You can now continue to improve the connectivity of the model by
%identifying metabolites that should be connected. A convenient way to get
%an overview of how connected the model is, and at the same time getting a
%lot of useful data, is to use gapReport. Note that it can take several
%hours to run, depending on the number of gaps in the model.
[noFluxRxns noFluxRxnsRelaxed subGraphs notProducedMets minToConnect...
    neededForProductionMat]=gapReport(newModel);

%You will see that 361/1140 reactions cannot carry flux. Remember that we
%allow for output of all metabolites, that's why it calculates 361 in both
%cases. 557/1195 metabolites cannot be synthesized from the precursors we have
%supplied. There are 7 subnetworks in the model, of which 1181/1195
%metabolites belong to the first one.
%
%It will also print something similar to:
%
%To enable net production of all metabolites, a total of 266 metabolites
%must be connected
%Top 10 metabolites to connect:
% 	1. Acyl-[acyl-carrier protein][s] (connects 65 metabolites)
% 	2. 1-Radyl-2-acyl-sn-glycero-3-phosphocholine[s] (connects 16 metabolites)
% 	3. G00146[s] (connects 13 metabolites)
% 	4. NAD+[s] (connects 12 metabolites)
% 	5. Tetrahydrofolate[s] (connects 11 metabolites)
% 	6. Methylselenic acid[s] (connects 10 metabolites)
% 	7. G00003[s] (connects 10 metabolites)
% 	8. G00009[s] (connects 8 metabolites)
% 	9. FAD[s] (connects 6 metabolites)
% 	10. Thiamin monophosphate[s] (connects 6 metabolites)

%This is a very useful way of directing the gap-filling tasks to where they
%are of the greatest use. In this case it says that in order to have net
%production of all metabolites in the model from the simple inputs that we have
%defined (557 metabolites can currently not have net production) we need to
%connect 226 unconnected metabolites in. However, by connecting only the
%top 10 in the list 157/557 (28%) of them will now be connected.
%If we look at the list we see that the top ranking one is
%Acyl-[acyl-carrier protein]. This is because acyl-carrier protein acts as
%a carrier for many fatty acids, but it is not synthesized in the model
%from its amino acid constituents. Therefore there cannot be net production
%of any of the metabolites involving it. One solution would be to add a
%reaction representing the synthesis of it (based on its amino acid
%composition). Another solution would be to add an uptake reaction for it.
%Metabolites 4,5,8, and 9 are all co-factors or involved in co-factor
%metabolism. It turns out that yeast cannot grow on only the substrates we
%have defined, but that it requires some other precursors for co-factor
%synthesis as well.

%Add uptake reactions for the minimal media constituents needed for yeast
%to grow.
[I J]=ismember({'4-Aminobenzoate';'Riboflavin';'Thiamine';'Biotin';'Folate';'Nicotinate';'Zymosterol';'Choline'},newModel.metNames);
[newModel addedRxns]=addExchangeRxns(newModel,'in',J);

%Rerun gapReport and use the output for targeting the gap-filling efforts.
%Note that only some info is printed; most of it is available in the output
%structures.
%Work like this in an iterative manner until the model is of sufficient
%quality.