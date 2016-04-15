%****RAVEN TOOLBOX MINI TUTORIAL***
%
% This small tutorial will use the Penicillium chrysogenum model to
% show some of the functionality of RAVEN Toolbox, with a focus on the 
% import/export functions and how to interpret the simulation data. 
% For a more detailed description of the individual functions, 
% see /doc/index.html.
%
% We'll be using a genome-scale model (GEM) of the filamentous fungi Penicillium
% chrysogenum. The model can be found in a Microsoft Excel file in raven.zip
% under the name iAL1006 v1.00.xlsx
%
% It is assumed that all files are in the current directory or in the Matlab
% path list.
%
% Rasmus Agren, 2013-08-06

%Import the model from Excel. This function performs a number of checks
%regarding the model structure (such as for incorrectly written equations or
%illegal characters). In this structure there is only one
%warning; that the formula for the metabolite LPE couldn't be parsed. The
%"false" flag imports a model with exchange reactions in their "closed"
%form. This makes the model unsuited for modelling, but it is useful for
%some quality control steps
model=importExcelModel('iAL1006 v1.00.xlsx',false)

%The Excel interface only works on Windows. On other systems you would need
%to import the model from SBML instead
%model=importModel('iAL1006 v1.00.xml',false)

%This function prints some properties of the model. The two "true" flags
%say that it should also list potential problems such as dead-end reactions
%or unconnected metabolites
printModelStats(model,true,true)

%As can be seen the model contains 1632 reactions, 1395 metabolites, and
%1006 genes

%Most modelling approaches using GEMs are based on mass balancing around
%the internal metabolites in the system. However, in order for the system
%to take up or excrete metabolites, some metabolites have been defined as
%unconstrained. In order to simulate something those metabolites have to be
%removed from the model. The function simplifyModel is a general-purpose
%function for making models smaller. This includes options such as
%grouping linear reactions and deleting reactions which cannot carry flux.
%Here we choose to delete the exchange metabolites, all reactions that are
%constrained to zero (mainly uptake of non-standard carbon sources), and
%all reactions that cannot carry flux (mainly reactions that were dependent
%on any of those non-standard carbons sources).
model=simplifyModel(model,true,false,true,true)

%As can be seen the model now contains only 1305 reactions, 1037 
%metabolites, and 1006 genes.

%Let's say that we want to calculate the theoretical yield of carbon
%dioxide from glucose as a first way of validating the model. The supplied
%model already allows for uptake of phosphate, sulfate, NH3, O2 and the co-factor
%precursors thiamin and pimelate.
%The setParam function is used for setting constraints, reversibility
%and objective function coefficients.

%Set the uptake of glucose to be no more than 1 mmol/gDW/h and no uptake of
%ethanol
model=setParam(model,'ub',{'glcIN' 'etohIN'},[1 0]);

%Set the objective for the simulation to maximize CO2 production
model=setParam(model,'obj',{'co2OUT'},1);

%We can then solve the problem using linear programming. The solveLP
%function takes a model and solves the linear programming problem defined
%by the constraints and the objective value coefficients.
sol=solveLP(model)

%If everything worked fine you should see a structure that contains the
%fields .f which is the negative of the objective value, .stat which is 1
%if the optimization terminated successfully and .x which contains the
%fluxes through each of the reactions in this particular solution of the
%problem. sol.x contains 1305 values which makes it rather difficult to
%interpret. A first approach is to look at only the reactions that
%transports metabolites in to or out from the system. The printFluxes
%function prints (parts of) the flux distribution to the console.

printFluxes(model, sol.x, true, 10^-7); %true means only print exchange fluxes

%We can see that the system took up one unit of glucose and 6 units of O2 while
%producing 6 units of CO2 and 6 units of H2O. Let's say that we want to
%understand more in detail how this happens. We can use the same function
%to plot all fluxes (above 10^-7)
printFluxes(model, sol.x, false, 10^-7);

%You will see that there are many reactions that have -1000 or 1000 flux.
%This is because there are loops in the solution. In order to clean up the
%solution we can minimize the sum of all the fluxes. This is done by setting 
%the third flag to solveLP to true (take a look at solveLP, there are other
%options as well)
sol=solveLP(model,1);
printFluxes(model, sol.x, false, 10^-7);

%You will see that there are now much fewer reactions that are active. Some
%of them will probably appear to be a little "weird" but that is because
%the ATP and NAD(P)H produced has to go somewhere in order to be mass balanced.

%If we want to study growth instead we simply change the objective to be
%production of biomass instead of CO2
model=setParam(model,'obj',{'bmOUT'},1);
sol=solveLP(model,1);
printFluxes(model, sol.x, true, 10^-7);

%We can see that the growth rate is 0.084/h and that the system now also
%requires sulfate,phosphate, NH3, thiamin and pimelate. Let's compare 
%this to the growth on ethanol instead of glucose. We use three times the 
%molar flux of ethanol since it contains 2 carbons rather than 6. 'eq' 
%means 'equal to' and sets the lower and upper bound to the same value.
modelETH=setParam(model,'eq',{'glcIN' 'etohIN'},[0 3]);
solETH=solveLP(modelETH,1);
printFluxes(modelETH, solETH.x, true, 10^-7);

%Say that we are interested in how metabolism changes between the two
%carbon sources. followChanged takes two flux distributions and lets you
%select which reactions to print. Here we show reactions that differ with
%more than 50%, has a flux higher than 0.5 mmol/gDW/h and an absolute
%difference higher than 0.5 mmol/gDW/h. 
followChanged(modelETH,sol.x,solETH.x, 50, 0.5, 0.5);

%There are 65 such reactions. By studying them you can start to get an idea
%about where the key changes occur. Visualization can help a lot in this
%regard. Say that we are particularly interested in how ATP metabolism
%changes. Then we can show only such reaction by writing,
followChanged(modelETH,sol.x,solETH.x, 30, 0.4, 0.4,{'ATP'});

%We can then see that on glucose ATP is generated in glycolysis but on
%ethanol it seems to have to do with acetate and so on. This allows you
%to look further and further down until you can understand the underlying
%flux redistributions that give rise to different phenotypes.

%The fluxes can also be visualized on a metabolic map. Green corresponds to
%reactions which are more used for growth on glucose, and red are reactions 
%which are more used for growth on ethanol. Open the GLCvsETH.pdf PDF to be
%able to zoom in on individual reactions.
load 'pcPathway.mat' pathway;
drawMap('Glucose vs ethanol',pathway,model,sol.x,solETH.x,modelETH,'GLCvsETH.pdf',10^-5);

