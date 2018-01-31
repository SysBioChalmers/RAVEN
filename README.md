# ![The RAVEN Toolbox 2](RAVEN2.png)
The RAVEN (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox 2 is a software suite for Matlab that allows for semi-automated reconstruction of genome-scale models (GEMs). It makes use of published models and/or KEGG, MetaCyc databases, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology.

If you are using RAVEN in any scientific work, please cite: [R. Agren, et. al, “The RAVEN Toolbox and Its Use for Generating a Genome-scale Metabolic Model for Penicillium chrysogenum,” PLoS Comput. Biol., vol. 9, no. 3, p. e1002980, Mar. 2013.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980).

> A manuscript describing RAVEN Toolbox 2 is currently being prepared. Citation details will therefore be updated in the near future.

Please report any technical issues and bugs [here](https://github.com/SysBioChalmers/RAVEN/issues). For other issues, please contact [Eduard Kerkhoven](https://github.com/edkerk).

-----
## Overview

[Releases](#releases)

[Installation](#installation)

[Tutorials](#tutorials)

[Hidden Markov Models (HMMs) for KEGG reconstruction](#hidden-markov-models-for-kegg-based-reconstruction)

[Development](#development-guidelines)

[Links](#links)

-----

## Releases
RAVEN can be installed via cloning the GitHub repository as per below or by downloading and extracting one of the a zipped [release](https://github.com/SysBioChalmers/RAVEN/releases). Please note that the releases do not always represent the most up to date version.

## Installation
### Required software
* A functional [MATLAB](mathworks.com/products/matlab.html) installation (version 2013b or later).

### Dependencies
* [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/5.15.0/stable/MATLAB%20interface/) (version 5.15 or higher), which is utilised for importing and exporting GEMs in SBML format. Note: not needed if [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) is installed.
* At least one solver for linear programming:
  * Preferred: [Gurobi Optimizer](http://www.gurobi.com/downloads/gurobi-optimizer) (version 7.5 or higher), academic license is available [here](https://user.gurobi.com/download/licenses/free-academic).
  * Alternative/legacy: [MOSEK](https://www.mosek.com/downloads/details/5/) (version 7 only), academic license is available [here](https://www.mosek.com/products/academic-licenses/).
  * If the user has [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) installed, it is possible to use the default COBRA solver (the one which is set by [_changeCobraSolver_](#cobra-toolbox).


### Instructions
#### RAVEN Toolbox
Obtain a RAVEN Toolbox in one of the following ways:
* In Terminal/Command Prompt, navigate to the desired installation directory and run the following Git command:
```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```
* Alternatively, download the latest [release](https://github.com/SysBioChalmers/RAVEN/releases) of RAVEN Toolbox as a ZIP file, and extracted to your favourite directory.

Once extracted, ensure that all other software dependencies (e.g. libSBML, Gurobi) are installed (see above for [list](#dependencies), below for [instructions](#libsbml). Then, open MATLAB and run the following command:
```matlab
cd('[location]/RAVEN/installation'))
checkInstallation
```
where ```[location]``` is the directory where you installed RAVEN.

This function checks the functionality for libSBML MATLAB API and solver software. It automatically recognises which solvers are installed and sets the first functional solver as the default RAVEN solver. The default RAVEN solver be changed any time by typing in Matlab:

```matlab
setRavenSolver('solverName')
```

Available solver names are ```gurobi```, ```mosek``` and ```cobra```.

In Unix-based systems _checkInstallation_ also checks the consistency of external binary programs. If these binaries are broken, they need to be re-compiled from their corresponding source codes. See the documentation for the corresponding software for more details.

#### libSBML
1. Download libSBML from the link [above](#dependencies) and install to your favourite directory.
2. In MATLAB, run the following command:

```matlab
addpath('[location]/libSBML-5.x.0-matlab')
savepath
```

where ```[location]``` is where you installed libSBML and ```5.x.0``` is your libSBML version (5.15.0 or higher).

#### Gurobi
1. Download from the link [above](#dependencies) and install Gurobi to your favourite location.
2. Make sure you obtained a [license](https://user.gurobi.com/download/licenses/free-academic) following instructions for [Windows](https://www.gurobi.com/documentation/7.5/quickstart_windows/retrieving_and_setting_up_.html), [Mac](https://www.gurobi.com/documentation/7.5/quickstart_mac/retrieving_and_setting_up_.html) or [Unix](https://www.gurobi.com/documentation/7.5/quickstart_linux/retrieving_and_setting_up_.html).
3. To install Gurobi in MATLAB, follow the instructions for [Windows](https://www.gurobi.com/documentation/7.5/quickstart_windows/matlab_setting_up_gurobi_f.html), [Mac](https://www.gurobi.com/documentation/7.5/quickstart_mac/matlab_setting_up_gurobi_f.html) or [Unix](https://www.gurobi.com/documentation/7.5/matlab_setting_up_gurobi_f.html).
4. Make sure that MATLAB remembers the Gurobi installation for next time, by running the following command: 

```matlab
savepath
```

#### Mosek
1. Download from the link [above](#dependencies) and install Mosek to your favourite location.
2. Make sure you obtained a [license](https://www.mosek.com/products/academic-licenses/) following [instructions](https://docs.mosek.com/8.0/install/installation.html#setting-up-the-license).
3. To install Mosek in MATLAB, follow [instructions](https://docs.mosek.com/8.0/toolbox/installation.html#id1). Note: the documentation mentions version 8, but RAVEN only works with version 7 of Mosek.
4.  Make sure that MATLAB remembers the Mosek installation for next time, by running the following command: 

```matlab
savepath
```

#### COBRA Toolbox
1. To gain access to functions from COBRA Toolbox, follow installation instructions provided [here](https://opencobra.github.io/cobratoolbox/latest/installation.html).
2. To use COBRA-specified solvers (e.g. open-source GLPK solver), configure COBRA and RAVEN with the following commands:

```matlab
changeCobraSolver('glpk')
setRavenSolver('cobra')
```

## Tutorials
Some tutorials highlighting basic RAVEN functionality can be found in the 'tutorial' folder in the installation directory.

## Hidden Markov Models for KEGG based reconstruction

Provided are pre-trained Hidden Markov Models (HMMs) for KEGG Orthology (KO) protein sets:

##### For RAVEN 2.0 and later

For _de novo_ reconstruction of a GEM, the RAVEN function _getKEGGModelForOrganism_ can use HMMs trained on KO protein sets. Provided are HMMs trained on KEGG Release 82.0. CD-HIT was used to obtain non-redundant representative KO protein sets thereby clustering proteins with the defined identity and overlap with the longest protein in the corresponding cluster threshold values. Multisequence alignment with MAFFT and training with HMMER 3.1b2 were then performed. The provided archives contain only pre-trained HMMs.

HMM sets can be downloaded **automatically** during GEM reconstruction from KEGG (set the *dataDir* parameter in *getKEGGModelForOrganism*). Alternatively, download links are provided below. The following HMM sets are available:

| dataDir | KEGG version | Phylogeny | Identity (%) | Overlap (%) |
|:-------:|:------------:|:---------:|:------------:|:-----------:|
|[euk100_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|eukaryota|100|90
|[euk90_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|eukaryota|90|90
|[euk50_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|eukaryota|50|90
|[prok100_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|prokaryota|100|90
|[prok90_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|prokaryota|90|90
|[prok50_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|82.0|prokaryota|50|90

##### For RAVEN 1.9.0 or older

HMMs were trained from KO protein sets, based on KEGG Release 58.1. Multisequence alignment was performed with ClustalW2, whereas HMMs were trained with HMMER 2.3. All the associated proteins were used in multisequence alignment and HMMs training. In addition to pre-trained HMMs, the archives also contain multisequence alignment data. The following HMM sets are available:

| Dataset | KEGG version | Phylogeny |
|:-------:|:------------:|:---------:|
|[eukaryota](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|58.1|eukaryota
|[prokaryota](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|58.1|eukaryota

## Development guidelines

Anybody is welcome to contribute to the development of RAVEN Toolbox, but please abide by the following guidelines.

### Urgent bugfix
* When fixing a bug in an existing function, make a separate branch from `master` and name the branch after the function you are fixing.
* Make commits to this branch while working on your bugfix. Note that bugfixes have to be backwards compatible.
* When you are certain that your bugfix works, make a pull request to the `master` branch. Also, see [Pull request](#pull-request) below.

### New features/functions
* For other development (not bugfixes, but for instance introducing new functions or new/updated features for existing functions): make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing.
* Make commits to this branch while developing. Aim for backwards compatibility, and try to avoid very new MATLAB functions when possible, to accommodate users with older MATLAB versions.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries, KEGG or MetaCyc database files, etc.
* `doc`: updating documentation (in `doc` folder) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: exportModel additional export to YAML
chore: update KEGG model to version 83.0
fix: optimizeProb parsing results from Gurobi
```
More detailed explanation or comments can be left in the commit description.

### Pull request
* No changes should be directly commited to the `master` or `devel` branches. Pull requests should be used.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* Typically, wait ~ 1 week before merging, to allow for other developers to inspect the pull request.
* A merge with the master branch typically invokes a new release (see [versioning](#versioning)).

### Versioning
RAVEN Toolbox follows [semantic versioning](https://semver.org/), and a `version.txt` file is updated with each release of the master branch.

## Links
For more systems biology related software and recently published genome-scale models from the Systems and Synthetic Biology group at Chalmers University of Technology, please visit the [Github page](https://github.com/SysBioChalmers). For more information and publications by the Systems and Synthetic Biology please visit [SysBio](www.sysbio.se).