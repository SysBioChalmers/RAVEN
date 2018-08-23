Changes after version 1.8.0 are tracked on GitHub, no detailed changelog is kept. Users are instead refered to: https://github.com/SysBioChalmers/RAVEN/compare/3d8d2223f796a70f033127e19c22b55db66ba155...master

Below is a non-exhaustive summary of changes made since version 1.8.0:
-Numerous bug-fixes related to many functions
-I/O support for SBML L3V1 FBCv2
-I/O support for Excel using Apache POI
-Support for writing COBRA(py) compatible YAML and TXT format
-Support for maintaining models in GitHub using exportForGit function
-Novel MetaCyc-based de novo reconstruction functionality
-Rename function names to prevent clashes with COBRA Toolbox
-ravenCobraWrapper to allow exchanging models between the two toolboxes
-getKEGGModelForOrganism is more efficient, implementing CD-HIT, MAFFT and HMMER3
-Flux scanning with enforced objective function (FSEOF) analysis included
-Include spontaneous reactions and transport reactions from MetaCyc

======================
The RAVEN Toolbox version 1.8.0
2016-04-15

Note: Since moving to GitHub a versioning system on the form vx.y.z has been adopted. E.g. version 1.08 corresponds to v.1.8.0

changelog v1.8.0
----------------------
This is the first release of RAVEN on Github with a change to the versioning system and corresponds to the earlier RAVEN v1.08. Installation instructions can be found in the install.txt file.

======================
-The RAVEN Toolbox version 1.08
2014-01-06

changelog v1.08
----------------------
-Fixed a bug in getEssentialRxns which caused it to misinterpret some problems as infeasible
-Fixed a bug in randomSampling where it could miss that a problem was infeasible
-Added a setting to convertToIrrev to only split some subset of the reversible reactions
-Fixed a bug in getModelFromHomology regarding quotation marks in the description field
-Fixed a bug in getKEGGModelForOrganism for some combination of settings
-Added the getRxnsInComp function
-Changed in checkRxn to hide an annoying warning
-Fixed a bug in copyToComps which caused the model to be renamed
-Added the cleanBadChars and cleanBadCharsInModel functions for ensuring that illegal characters are dealt with before generation of SBML files
-Fixed a bug in checkRxns which caused an error for some combination of inputs
-Fixed a bug in getBlast when using multiple reference files
-Fixed a bug in getModelFromHomology for some combination of parameters
-Added a setting to constructS for printing more detailed warnings/errors. The purpose was to be able to use it as the general parsing function for importExcelModel
-Added the checkModelStruct function
-Fixed a bug in predictLocalization for some models without gene complexes
-Fixed a bug in removeBadRxns
-Changed in parseFormulas so that all elements can be used in metabolite composition. This also required changes in getElementalBalance, printFluxes, printModelStats and removeBadRxns
-Fixed a bug in predictLocalization which could result in an Index out of bounds error in some cases when using only two compartments
-Added an option to use different transportCosts for each metabolite to predictLocalization
-Rewrote how import/output of models work. Instead of using SBMLFromExcel followed by importModel, it is now possible to use importExcelModel directly. Most of the checks which were performed in SBMLFromExcel have now been moved to checkModelStruct, which can then also be used for any model structure. checkModelStruct is also called from importModel/importExcelModel. The report about dead-end reactions and similar issues which do not relate directly to the model structure as such have been moved to printModelStats. The exportModel function has been rewritten and is now equally exact as first exporting to Excel format and then to SBML from there. The overall ambition with the rewrite was to modularize the input/output and make it easier to perform quality controls of model structures without going via Excel. SBMLFromExcel is now just a wrapper for importExcelModel+exportModel. Another reason was to get rid of largely redundant code in SBMLFromExcel/importModel/exportModel. This work resulted in some small changes to the Excel model format and the Matlab model structure.
	-The "SBO TERM" column for reactions is no longer used and no SBO terms are incorporated in the SBML file
	-Differences in upper/lower case when matching metabolites and reactions are now considered errors rather than warnings
	-The model structure has an added "annotation" field used for meta data
	-Several Miriam terms can now be delimited with ";"
	-"Gene ID 1", "Gene ID 2" and "KEGG MAPS" in the  GENE sheet have been removed and "MIRIAM" has been added
	-"GO TERM" has been removed from the COMPS sheet and "MIRIAM" has been added
	-"MIRIAM" had been added to the RXNS sheet
	-The model structure can now contain the fields "rxnComps", "rxnMiriams", "compMiriams", "geneMiriams" and "geneComps"
	-The model field geneShortName is no longer used
	-The naming of columns in Excel sheets have changed to be "ID" & "NAME" for RXNS and METS, "ID" & "DESCRIPTION" for
	 MODEL, "NAME" for GENES and "ABBREVIATION" & "NAME" for COMPS. The old captions still work as well.
-Added support for importing MS Excel files on non-Windows operating systems. This makes use of Oracle POI
-Changed in getKEGGModelForOrganism so that a model for an organism which is already in KEGG could be generated without any sequence data
-Changed in getAllSubGraphs so that it doesn't rely on the Graph toolbox
-Fixed a bug in guessComposition which gave errors for some polymers
-Fixed a bug in addRxns


======================
-The RAVEN Toolbox version 1.7
2013-07-16

changelog v1.07
----------------------
-Fixed a bug in runINIT when using metabolomics data without free excretion
-Added the checkRxn function
-Fixed a bug in drawMap regarding non-mapped reactions
-Fixed a bug in drawMap regarding the cutOff value
-Fixed in checkSolution so that it throws an error if the licence file is missing
-Changed in exportToExcelFormat so that it doesn't create the empty default sheets "Sheet1", "Sheet2", and "Sheet3"
-Fixed a bug in removeMets which caused it do not delete compMiriams if removeComps=true
-Fixed a bug in fillGaps which caused the function to fail to throw an error if the constrains were infeasible also when using the full template model
-Changed in scoreModel so that it throws an error if arrayData only contains one tissue/celltype
-Changed in solveLP, getMILPParams, and runINIT to make RAVEN work with Mosek 7

======================
-The RAVEN Toolbox version 1.06
2013-05-20

changelog v1.06
----------------------
-Added the option to print a report in removeBadRxns
-Fixed a bug in setParam when using logical vectors as input
-Added the option to ignore internal contraints in makeSomething/consumeSomething/removeBadRxns
-Added a check to getINITModel that the model should be in the closed form
-Fixed a bug in SBMLFromExcel which led to that InChI strings were ignored if no Miriam identifier was used
-The installation script (install.m) was removed. All it did was to add the paths to Mosek and libSBML to the the Matlab path list. I felt that it was better that the user just follows the installation directions for each software instead. Instead the function checkInstallation was added, which tests if everything is installed
-Fixed a bug in runINIT when presentMets is a row vector instead of a column vector

======================
-The RAVEN Toolbox version 1.05
2013-05-16

changelog v1.05
----------------------
-Changed so that removeBadRxns also looks at consumption of metabolites
-Changed in solveLP so that it prints a warning if the minimization of the sum of fluxes fails. It then returns the basic solution, but with exit flag -2. This required a small change in getEssentialRxns as well
-Fixed a bug in runINIT regarding compartment naming
-Rewrote solveLP for clarity and speed. No change in functionality, but more numerically stable
-Added the checkSolution function for checking that a solution is feasible/optimal. These checks in solveLP, getMinNrFluxex and runINIT were removed, and those functions call on checkSolution instead
-Fixed a bug in getINITModel where it crashed if no taskStructure was supplied

======================
-The RAVEN Toolbox version 1.04
2013-05-11

changelog v1.04
----------------------
-Fixed a problem with the iAL1006 model in the tutorial, which resulted in that it didn't fit with the code in tutorial0.m. The reason was that all input reactions were open by default, which should not have been the case
-Fixed an issue in the install script regarding the new version of libSBML

======================
-The RAVEN Toolbox version 1.03
2013-04-22

changelog v1.03
----------------------
-Added a check to SBMLFromExcel that the Excel file exists
-Fixed a bug in fitTasks which resulted in an error if the same id was used for both model and reference model
-Changed so that getMinNrFluxes only accepts negative reaction scores, as there was a problem with the implementation otherwise
-Added a check to fitTasks that all scores are negative
-Changed printModelStats to be faster. This came with a change in how unique reactions are defined. Before, a transport reaction like "A[c] (<)=> A[m]" was not included in the count because the metabolites cancelled in a non-compartmentalized model. Now it will count, but only once per transported metabolite ("A[c] (<)=> A[p]" will not give rise to another unique reaction)
-Added a setting to printModelStats to print a summary of potential issues with a model (dead-end reactions/metabolites, balancing problems and some other stuff). Can be used as a faster, but much less comprehensive, version of gapReport
-Changed so that getINITModel can accept a params structure for the fitTasks step, and not only for the INIT algorithm
-Changed so that getINITModel identifies dead-end reactions only based on topology and not by actually testing if flux through the reaction is possible. This was done because the original implementation takes several hours longer and the difference tends to be very small
-Changed so that the model returned from getINITModel contains the exchange reactions from the reference model. Before they tended to be deleted since they were most often dead ends in the closed model, and because they had negative scores (if not associated with any genes). Exchange reactions which involve metabolites that are not in the output model are still excluded
-Fixed a bug in randomSampling which could wrongly cause the model to be regarded as infeasible
-Added the addGenes function
-Fixed a bug in importModel regarding parsing of gene miriams codes
-Changed in solveLP so that it can make use of hot-start
-Changed in haveFlux so that it can make use of hot-start
-Changed in getAllowedBounds so that it can make use of hot-start
-Changed in checkTasks so that it returns the essential reactions for each task, rather than for the whole set of tasks
-changed in getEssentialRxns so that it can make use of hot-start
-Changed in fitTasks so that it returns the included reactions for each task, rather than for the whole set of tasks
-Changed in getINITModel so that it alse returns the essential reactions and the reactions included for gap-filling for each task, and not only for the whole set of tasks

======================
-The RAVEN Toolbox version 1.02
2013-03-05

changelog v1.02
----------------------
-Added a check for feasibility to solveQP
-Added a models simplification step to randomSampling for speed reasons
-Added an option in getElementalBalance to print warnings if reactions cannot be parsed
-Fixed a bug in getElementalBalance which caused it to only print some of the unbalanced reactions if that setting was chosen
-Added the option to printFluxes to print some details regarding elemental balancing
-Changed in getIndexes to accept character arrays and not only cell arrays
-Changed in printFluxes to make it faster
-Changed in randomSampling so that it doesn't rely on the Matlab Statistical Toolbox 
-Changed in randomSampling to speed up the identification of reactions involved in loops
-printModel was rewritten to be a wrapper function for printFluxes, since the functionality was largely redundant
-Fixed a bug in exportToExcelFormat regarding Inf bounds
-Changed in getElementalBalance so that it returns the left/right hand side even if the reaction as a whole cannot be parsed
-The default cutoff in haveFlux was lowered to 1e-8 to reduce the risk of removing functional reactions
-Changed in setParam so that it can accept indexes and not only reaction ids
-Fixed a bug in install.m
-Changed in checkTasks so that it allows new metabolites in added reactions
-Fixed a bug in fitTasks relating to the new way to store compartment for metabolites

======================
-The RAVEN Toolbox version 1.01
2013-02-07

changelog v1.01
----------------------
-Added some checks to getEssentialRxns to make it less sensitive to numerical issues
-Fixed a bug in mergeModels regarding the compOutside field
-Changed in getModelFromHomology so that the output model has matching "id" and "description" fields
-Changed in SBMLFromExcel so that reactions can be written with only one side defined ("<=> met1" or "met1 => ")
-Fixed a bug in exportToExcelFormat regarding the compOutside field
-Fixed a bug in consumeSomething/produceSomething when the only possible solution was to consume/produce several metabolites at the same time
-Added the copyToComps function. Is a utility function for copying reactions to other compartments
-Added the getMetsInComp function. Is a utility function for getting which metabolites are in a given compartment
-Changed so that not all models need to have an unconstrained field in mergeModels
-Change the way compartments are stored in the structure. "comps" now represent an abbreviation of the compartment. "compNames" are the names of the compartments. "compOutside" is the abbreviation of the outside compartments. "metComps" is the index in comps. Since "comps" now can be more descriptive it is used rather than "compNames" when displaying compartment information. This change affects many of the RAVEN functions and it is not totally backward compatible