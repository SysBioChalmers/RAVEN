% tutorial0
%	This is a short introduction that shows how to load a genome-scale
%	metabolic model (GEM), set reaction constraints, objective function and
%	perform an optimization through flux balance analysis (FBA). The
%	resulting fluxes are visualized and exported to a PDF file.
%   For a more detailed description of the individual functions, see
%   [raven_directory]/doc/index.html.
%   A GEM for the filamentous fungus Penicillium chrysogenum is used in
%   this tutorial. The model can be found in a Microsoft Excel file under
%   the name "iAL1006 v1.00.xlsx" and in SBML file "iAL1006 v1.00.xml".
%   See Exercise 0 in "RAVEN mini tutorial.docx" for more details.
%
%	Simonas Marcisauskas, 2019-10-21
%

%Import the model from Excel. This function performs a number of checks
%regarding the model structure (such as for incorrectly written equations
%or illegal characters). In this structure there is only one warning; that
%the formula for the metabolite LPE could not be parsed. The "false" flag
%imports a model with exchange reactions in their "closed" form. This makes
%the model unsuited for modelling, but it is useful for some quality
%control steps.
model=importExcelModel('iAL1006 v1.00.xlsx',false);

%The Excel interface is supposed to work in all the systems (Windows, Unix,
%macOS), but upon any problems, the model can be imported from SBML format
%instead. However, in such case the user would not be able to run
%tutorial1, tutorial2 and tutorial3, since these tutorials involve the
%editing of Excel files. Run the command below (remove "%" sign) instead,
%if having such problem:
%model=importModel('iAL1006 v1.00.xml',false);

%The following function prints some properties of the model. The two "true"
%flags say that it should also list potential problems such as dead-end
%reactions or unconnected metabolites.
printModelStats(model,true,true);

%As can be seen the model contains 1632 reactions, 1395 metabolites and
%1006 genes

%Most modelling approaches using GEMs are based on the mass balancing
%around the internal metabolites in the system. However, in order for the
%system to uptake or excrete metabolites, some metabolites have been
%defined as "unconstrained". In order to simulate something, those
%metabolites have to be removed from the model. The function simplifyModel
%is a general-purpose function for making models smaller. This includes the
%options such as grouping linear reactions and deleting reactions which
%cannot carry flux. Here it is chosen to delete the exchange metabolites,
%all reactions that are constrained to zero (mainly uptake of non-standard
%carbon sources), and all reactions that cannot carry flux (mainly
%reactions that were dependent on any of those non-standard carbons
%sources).
model=simplifyModel(model,true,false,true,true);

%As can be seen the model now contains only 1305 reactions, 1037
%metabolites and 1006 genes

%As a first way of validating the model, calculate the theoretical yield of
%carbon dioxide from glucose. The supplied model already allows for uptake
%of phosphate, sulfate, NH3, O2 and the co-factor precursors thiamin and
%pimelate. The setParam function is used for setting constraints,
%reversibility and objective function coefficients.
%Set the uptake of glucose to be no more than 1 mmol/gDW/h and no uptake of
%ethanol.
model=setParam(model,'ub',{'glcIN' 'etohIN'},[1 0]);

%Set the objective for the simulation to maximize CO2 production
model=setParam(model,'obj',{'co2OUT'},1);

%The problem can now be solved using linear programming. The solveLP
%function takes a model and solves the linear programming problem defined
%by the constraints and the objective value coefficients.
sol=solveLP(model);

%If everything worked fine one should see a structure that contains the
%fields .f which is the negative of the objective value, .stat which is 1
%if the optimization terminated successfully and .x which contains the
%fluxes through each of the reactions in this particular solution of the
%problem. sol.x contains 1305 values which makes it rather difficult to
%interpret. A first approach is to look at only the reactions that
%transports metabolites in to or out from the system. The printFluxes
%function prints (parts of) the flux distribution to the console.
printFluxes(model, sol.x, true, 10^-7); %true means only print exchange fluxes

%One can see that the system took up one unit of glucose and 6 units of O2
%while producing 6 units of CO2 and 6 units of H2O. To get more detailed
%insight how this happens, one can use the same function to plot all fluxes
%(above 10^-7).
printFluxes(model, sol.x, false, 10^-7);

%The results show many reactions that have -1000 or 1000 flux. This is
%because there are loops in the solution. In order to clean up the solution
%one can minimize the sum of all the fluxes. This is done by setting the
%third flag to solveLP to true (take a look at solveLP, there are other
%options as well).
sol=solveLP(model,1);
printFluxes(model, sol.x, false, 10^-7);

%Now there are much fewer reactions that are active. Some of them will
%probably appear to be a little "weird" but that is because the ATP and
%NAD(P)H produced has to go somewhere in order to be mass balanced.

%If one wants to study the growth instead, change the objective to be the
%production of biomass instead of CO2
model=setParam(model,'obj',{'bmOUT'},1);
sol=solveLP(model,1);
printFluxes(model, sol.x, true, 10^-7);

%The results show that the growth rate is 0.084/h and that the system now
%also requires sulfate, phosphate, NH3, thiamin and pimelate. Compare this
%to the growth on ethanol instead of glucose. Use three times the molar
%flux of ethanol since it contains 2 carbons rather than 6. 'eq' means
%'equal to' and sets the lower and upper bound to the same value.
modelETH=setParam(model,'eq',{'glcIN' 'etohIN'},[0 3]);
solETH=solveLP(modelETH,1);
printFluxes(modelETH, solETH.x, true, 10^-7);

%Investigate how metabolism changes between the two carbon sources.
%followChanged takes two flux distributions and lets the user select which
%reactions should be printed. Here the reactions are shown that differ with
%more than 50%, has a flux higher than 0.5 mmol/gDW/h and an absolute
%difference higher than 0.5 mmol/gDW/h.
followChanged(modelETH,sol.x,solETH.x, 50, 0.5, 0.5);

%There are 65 such reactions. By studying them one can start to get an idea
%about where the key changes occur. Visualization can help a lot in this
%regard. One can investigate how ATP metabolism changes by running the
%following command:
followChanged(modelETH,sol.x,solETH.x, 30, 0.4, 0.4,{'ATP'});

%See that on glucose ATP is generated in glycolysis but on ethanol it seems
%to have to do with acetate and so on. This allows the user to look further
%and further down until one understands the underlying flux redistributions
%that give rise to different phenotypes.

%The fluxes can also be visualized on a metabolic map. Green corresponds to
%reactions which are more used for growth on glucose, and red are reactions
%which are more used for growth on ethanol. Open the "GLCvsETH.pdf" file to
%be able to zoom in on individual reactions.
load 'pcPathway.mat' pathway;
drawMap('Glucose vs ethanol',pathway,model,sol.x,solETH.x,modelETH,'GLCvsETH.pdf',10^-5);
