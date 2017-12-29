# The RAVEN Toolbox
The RAVEN (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite for Matlab that allows for semi-automated reconstruction of genome-scale models (GEMs). It makes use of published models and/or KEGG, MetaCyc databases, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology.

If you are using RAVEN in any scientific work, please cite: [R. Agren, et. al, “The RAVEN Toolbox and Its Use for Generating a Genome-scale Metabolic Model for Penicillium chrysogenum,” PLoS Comput. Biol., vol. 9, no. 3, p. e1002980, Mar. 2013.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980).

> A manuscript describing RAVEN Toolbox 2.0 is currently being prepared. Citation details will therefore be updated in the near future.

Please report any technical issues and bugs [here](https://github.com/SysBioChalmers/RAVEN/issues). For other issues, please contact [Eduard Kerkhoven](https://github.com/edkerk).

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
  * If the user has [COBRA Toolbox](https://github.com/opencobra/cobratoolbox) installed, it is possible to use the default COBRA solver (the one which is set by _changeCobraSolver_).


### Instructions
#### RAVEN Toolbox
Obtain a RAVEN Toolbox in one of the following ways:
* In Terminal/Command Prompt, navigate to the desired installation directory and run the following Git command:
```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```
* Alternatively, download the latest [release](https://github.com/SysBioChalmers/RAVEN/releases) of RAVEN Toolbox as a ZIP file, and extracted to your favourite directory.

Once extracted, ensure that all other software dependencies (e.g. libSBML, Gurobi) are installed (see above for [list](#dependencies), below for [instructions](#libSBML). Then, open MATLAB and run the following command:
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
2. Make sure you obtained a [license](https://www.mosek.com/products/academic-licenses/) following [instructions] (https://docs.mosek.com/8.0/install/installation.html#setting-up-the-license).
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

## Pre-trained Hidden Markov Models (HMMs) for KEGG Orthology (KO) protein sets
_For RAVEN 2.0_

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

_For RAVEN 1.9.0 or older_

HMMs were trained from KO protein sets, based on KEGG Release 58.1. Multisequence alignment was performed with ClustalW2, whereas HMMs were trained with HMMER 2.3. All the associated proteins were used in multisequence alignment and HMMs training. In addition to pre-trained HMMs, the archives also contain multisequence alignment data. The following HMM sets are available:

| Dataset | KEGG version | Phylogeny |
|:-------:|:------------:|:---------:|
|[eukaryota](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|58.1|eukaryota
|[prokaryota](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip)|58.1|eukaryota


## Links
For more systems biology related software and recently published genome-scale models from the Systems and Synthetic Biology group at Chalmers University of Technology, please visit the [Github page](https://github.com/SysBioChalmers). For more information and publications by the Systems and Synthetic Biology please visit [SysBio](www.sysbio.se).