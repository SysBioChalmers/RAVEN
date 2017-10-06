# The RAVEN Toolbox
The RAVEN (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite for Matlab that allows for semi-automated reconstruction of genome-scale models (GEMs). It makes use of published models and/or KEGG, MetaCyc databases, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology

If you are using RAVEN in any scientific work, please cite: [R. Agren, et. al, “The RAVEN Toolbox and Its Use for Generating a Genome-scale Metabolic Model for Penicillium chrysogenum,” PLoS Comput. Biol., vol. 9, no. 3, p. e1002980, Mar. 2013.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980).

Please report any technical issues and bugs [here](https://github.com/SysBioChalmers/RAVEN/issues). For other issues, please contact [Eduard Kerkhoven](https://github.com/edkerk).

## Releases
RAVEN can be installed via cloning the GitHub repository as per below or by downloading and extracting one of the a zipped [release](https://github.com/SysBioChalmers/RAVEN/releases). Please note that the releases do not always represent the most up to date version.

## Installation
1. In Terminal/Command Prompt, navigate to the desired installation directory and run the following Git command:

```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```

Alternatively, the whole RAVEN repository can be downloaded as a ZIP file and extracted in the desired installation directory.

Once the RAVEN repository is downloaded, add its directory with subdirectories to the Matlab path. To do so, open Matlab and type:

```matlab
addpath(genpath('path/to/RAVEN'))
```

'path/to/' is the absolute path to the RAVEN Toolbox installation directory

2. Install [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/5.15.0/stable/MATLAB%20interface/) (v. 5.15 or higher), which is utilised for importing and exporting GEMs in SBML format. Install it and add to the Matlab path:

```matlab
addpath(genpath('path/to/libSBML'))
```

The installation of libSBML MATLAB API is not needed if the user has installed [the COBRA Toolbox](https://github.com/opencobra/cobratoolbox/)

3. Install solver software. For simulation and gap filling purposes at least one solver must be installed. The following solvers are compatible with RAVEN:
- [Gurobi Optimizer](http://www.gurobi.com/downloads/gurobi-optimizer). Academic licenses are available [here](https://user.gurobi.com/download/licenses/free-academic). It is strongly recommended to use this software as we support up-to-date versions of it. Install the solver, activate the license and add the installation directory to the Matlab path:

```matlab
addpath(genpath('path/to/gurobi'))
```

- [MOSEK](https://www.mosek.com/downloads/). Academic licenses are available [here](https://www.mosek.com/products/academic-licenses/). Only version 7 is compatible with the RAVEN Toolbox. Install the solver, activate the license and add the installation directory to the Matlab path:

```matlab
addpath(genpath('path/to/mosek'))
```

- If the user has the COBRA Toolbox installed, it is possible to use the default COBRA solver (the one which is set by _changeCobraSolver_). In such case the user only needs to ensure that the solver is already set by _changeCobraSolver_.

4. In Matlab type:

```matlab
checkInstallation
```

This function checks the functionality for libSBML MATLAB API and solver software. It automatically recognises which solvers are installed and sets the first functional solver as the default RAVEN solver. The default RAVEN solver be changed any time by typing in Matlab:

```matlab
setRavenSolver('solverName')
```

Available solver names are 'gurobi', 'mosek' and 'cobra'

In Unix-based systems _checkInstallation_ also checks the consistency of external binary programs. If these binaries are broken, they need to be re-compiled from their corresponding source codes. See the documentation for the corresponding software for more details


To have Matlab remember the above installations please save the current path to ’pathdef.m’ *in your MATLAB startup directory*. For instance on Mac OS X:

```matlab
savepath '~/Documents/MATLAB/pathdef.m'
```

*It is highly recommended to set set startup directory of Matlab to a directory where you have write access (Preferences->General->start up directory).*

## Tutorials
Some tutorials highlighting basic RAVEN functionality can be found in the 'tutorial' folder in the installation directory.

## Pre-trained Hidden Markov Models (HMMs) for KEGG Orthology (KO) protein sets
_For RAVEN > 1.9.0_

HMMs were trained from KO protein sets, based on KEGG Release 82.0. CD-HIT was used to obtain non-redundant representative KO protein sets thereby clustering proteins with the defined identity and overlap with the longest protein in the corresponding cluster threshold values. Multisequence alignment with MAFFT and training with HMMER 3.1b2 were then performed. The archives contain only pre-trained HMMs. Such HMM sets can be downloaded automatically during genome-scale metabolic model reconstruction from KEGG (see *dataDir* parameter in *getKEGGModelForOrganism*). The download links for HMM sets are also included below and in BioMet ToolBox. The following HMM sets are available:	
- [euk100_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk100_kegg82.zip). CD-HIT was used for all Eukaryotic proteins and the following threshold values: identity 100 %, overlap 90 %.

- [euk90_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk90_kegg82.zip). CD-HIT was used for euk100_kegg82 dataset and the following threshold values: identity 90 %, overlap 90 %.

- [euk50_kegg82](http://biomet-toolbox.org/tools/downloadable/files/euk50_kegg82.zip). CD-HIT was used for euk90_kegg82 dataset and the following threshold values: identity 50 %, overlap 90 %.

- [prok100_kegg82](http://biomet-toolbox.org/tools/downloadable/files/prok100_kegg82.zip). CD-HIT was used for all Prokaryotic proteins and the following threshold values: identity 100 %, overlap 90 %.

- [prok90_kegg82](http://biomet-toolbox.org/tools/downloadable/files/prok90_kegg82.zip). CD-HIT was used for prok100_kegg82 dataset and the following threshold values: identity 90 %, overlap 90 %.

- [prok50_kegg82](http://biomet-toolbox.org/tools/downloadable/files/prok50_kegg82.zip). CD-HIT was used for prok90_kegg82 dataset and the following threshold values: identity 50 %, overlap 90 %.


_For RAVEN 1.9.0 or older_

HMMs were trained from KO protein sets, based on KEGG Release 58.1. Multisequence alignment was performed with ClustalW2, whereas HMMs were trained with HMMER 2.3. All the associated proteins were used in multisequence alignment and HMMs training. In addition to pre-trained HMMs, the archives also contain multisequence alignment data. The following HMM sets are available:

- [eukaryota](http://biomet-toolbox.org/tools/downloadable/files/eukaryota.zip). Contains HMMs, trained from Eukaryotic organisms proteome.

- [prokaryota](http://biomet-toolbox.org/tools/downloadable/files/prokaryota.zip). Contains HMMs, trained from Prokaryotic organisms proteome.


## Links
For more information on software connected to Genome Scale models please visit the Systems Biology [Github page](https://github.com/SysBioChalmers). For information and publications by the Systems Biology department please visit [SysBio](www.sysbio.se).

## License (MIT)
Copyright (c) 2016 SysBio Chalmers

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.