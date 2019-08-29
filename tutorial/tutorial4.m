%This exercise is about creating a model from KEGG, based on protein
%sequences in a FASTA file, and doing some functionality checks on the
%model. The example case is for the yeast Saccharomyces cerevisiae. This
%tutorial is more of a showcase than the other three, and it's main purpose
%is to serve as a scaffold if you would like to reconstruct a GEM for your
%own organism.
%
%Simonas Marcisauskas, 2019-08-28
%

%Start by downloading trained Hidden Markov Models for eukaryotes. This can
%be done automatically or manually from BioMet ToolBox webpage. In this
%tutorial we will consider automatic download by picking euk100_kegg82
%archive. See the documentation in GitHub for more information regarding
%preparation of such archive.

%This creates a model for S. cerevisiae. The parameters are set to exclude
%general or unclear reactions and reactions with undefined stoichiometry.
%Type "help getKEGGModelForOrganism" to see what the different parameters
%are for. This process takes up to 20-30 minutes in macOS, Unix systems and
%40-50 minutes in Windows, depending on your hardware and the size of
%target organism proteome
model=getKEGGModelForOrganism('sce','sce.fa','euk100_kegg87','output',false,false,false,false,10^-30,0.8,0.3,-1);

%As you can see the resulting model contains (around) 1501 reactions, 1513
%metabolites and 813 genes. Small variations are possible since it is an
%heuristic algorithm.
model

%A first control is that your model shouldn't be able to produce any
%metabolites without uptake of some metabolites. This commonly happens when
%metabolites have a different meaning in different reactions. The best way
%to find such reactions is to run makeSomething and analyze the resulting
%solution for bad reactions. An automated approach is to use removeBadRxns,
%which tries to do the same thing in an automated manner. Please type "help
%removeBadRxns" for details.
[newModel, removedRxns]=removeBadRxns(model);

%You will see an error about that H+ can be made even if no reactions were
%unbalanced. Protons are particularly problematic since it's rather
%arbitary at which pH the formulas are written for. For the purpose of this
%analysis we can ignore protons and try to fix it later.
[newModel, removedRxns]=removeBadRxns(model,1,{'H+'},true);

%Only one reaction was removed because it enabled the model to produce
%something from nothing. Since its only one reaction, it might be
%worthwhile to look into this in more detail.
removedRxns

%If you look up this reaction it up in KEGG you will find that  it is a
%general polymer reaction. You might want to look at the flux distributions
%in more detail to try to find out if there is any better alternative to
%delete. Use makeSomething to do this
[fluxes, metabolite]=makeSomething(model,{'H+'},true);
model.metNames(metabolite)

%The model could produce H2O using the following reactions
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n')

%That resulted in quite a lot of fluxes to look through. Maybe it's easier
%if we exclude the ones that were elementally balanced.
balanceStructure=getElementalBalance(model);

%We know that hydrogen balancing is a bit tricky, so let's only look at the
%ones unbalanced for oxygen (since water was produced)
goodOnes=balanceStructure.leftComp(:,6)==balanceStructure.rightComp(:,6);
printFluxes(removeReactions(model,goodOnes), fluxes(~goodOnes), false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n')

%We still got a good number of reactions. Let's leave only the reactions
%which involve amylose or starch, from one of the problematic reactions
%that we identified earlier.
printFluxes(model, fluxes, false, [], [],'%rxnID (%rxnName):\n\t%eqn: %flux\n',{'Amylose';'Starch'});

%We got two elementally unbalanced reactions, including, the reaction
%which was identified by removeBadRxns. When looking to these reactions
%closer, one can notice the contradiction between the reactions. The first
%one shows that amylose and starch are interconvertible, while the second
%reaction shows that amylose contains one less glucose unit than starch.
%This type of general reactions are problematic and should be fixed
%manually. We therefore choose to trust removeBadRxns and delete R02110.
model=removeReactions(model,'R02110');

%The model can no longer make something from nothing. Can it consume
%something without any output?
[solution, metabolite]=consumeSomething(model,{'H+'},true);
model.metNames(metabolite)

%Nope, so that was good. Let's add some uptakes and see what it can
%produce.
[I, J]=ismember({'D-Glucose';'H2O';'Orthophosphate';'Oxygen';'NH3';'Sulfate'},model.metNames);
[model, addedRxns]=addExchangeRxns(model,'in',J);

%Check which metabolites can be produced given these uptakes. The
%canProduce function allows for output of all metabolites. This won't
%happen in the real cell, but it is very useful for functionality testing
%of the model. Once it is functional you can add excretion reactions based
%on evidence.
I=canProduce(model);

sum(I)/numel(model.mets)
%You can see that 27% of the metabolites could be synthesized. It is not
%directly clear whether this is a high or low number, many metabolites
%should not be possible to synthesize from those simple precursors.

%We can try to fill gaps using the full KEGG model to see if that gives a
%significantly higher number
keggModel=getModelFromKEGG([],false,false,false,false);

%The KEGG model is associated to more than 3,000,000 genes. They won't be
%used for the gapfilling, so we remove them to make this a little faster
keggModel=rmfield(keggModel,'genes');
keggModel=rmfield(keggModel,'rxnGeneMat');

%We've already seen that there are some unbalanced reactions in KEGG. Let's
%only use the balanced ones for the gap filling
balanceStructure=getElementalBalance(keggModel);
keggModel=removeReactions(keggModel,balanceStructure.balanceStatus~=1,true,true);

%fillGaps with these settings will try to include reactions in order to
%have flux through all reactions in the model. There are other settings as
%well. The first flag says that production of all metabolites should be
%allowed.
params.relGap=0.6; %Lower number for a more exhaustive search
params.printReport=true;
[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(model,keggModel,true,false,false,[],params);

%We see that we could connect 80 reactions (newConnected) by including 69
%reactions from the KEGG model (addedRxns). Those should of course be
%checked manually to see that they exist in yeast, but let's say that we
%have done so now.

%You can now continue to improve the connectivity of the model by
%identifying metabolites that should be connected. A convenient way to get
%an overview of how connected the model is, and at the same time getting a
%lot of useful data, is to use gapReport. Note that it can take several to
%many hours to run, depending on the number of gaps in the model.
[noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat]=gapReport(newModel);

%You will see that 486/1574 reactions cannot carry flux. Remember that we
%allow for output of all metabolites, that's why it calculates 486 in both
%cases. 687/1527 metabolites cannot be synthesized from the precursors we
%have supplied. There are 8 subnetworks in the model, of which 1511/1527
%metabolites belong to the first one.
%
%It will also print something similar to:
%
%To enable net production of all metabolites, a total of 304 metabolites
%must be connected
%Top 10 metabolites to connect:
% 	1. 2,3-Dehydroacyl-CoA[s] (connects 144 metabolites)
% 	2. 7,8-Dihydroneopterin[s] (connects 15 metabolites)
% 	3. Thiamin diphosphate[s] (connects 13 metabolites)
% 	4. Ergosterol[s] (connects 13 metabolites)
% 	5. Methylselenic acid[s] (connects 12 metabolites)
% 	6. G00003[s] (connects 12 metabolites)
% 	7. Androstenediol[s] (connects 11 metabolites)
% 	8. G00146[s] (connects 11 metabolites)
% 	9. 1-Alkyl-2-acylglycerol[s] (connects 10 metabolites)
% 	10. Progesterone[s] (connects 9 metabolites)

%This is a very useful way of directing the gap-filling tasks to where they
%are of the greatest use. In this case it says that in order to have net
%production of all metabolites in the model from the simple inputs that we
%have defined (687 metabolites can currently not have net production) we
%need to connect 304 unconnected metabolites in. However, by connecting
%only the top 10 in the list 209/687 (30%) of them will now be connected.
%If we look at the list we see that the top ranking one is
%2,3-Dehydroacyl-CoA, which is involved in the elongation of many different
%fatty acid chains. Metabolites 2 and 3 are co-factors or involved in co-
%factor metabolism. It turns out that yeast cannot grow on only the
%substrates we have defined, but that it requires some other precursors for
%co-factor synthesis as well.

%Add uptake reactions for the minimal media constituents needed for yeast
%to grow.
[I, J]=ismember({'4-Aminobenzoate';'Riboflavin';'Thiamine';'Biotin';'Folate';'Nicotinate';'Zymosterol';'Choline'},newModel.metNames);
[newModel, addedRxns]=addExchangeRxns(newModel,'in',J);

%Rerun gapReport and use the output for targeting the gap-filling efforts.
%Note that only some info is printed; most of it is available in the output
%structures. Work like this in an iterative manner until the model is of
%sufficient quality.
