# RAVEN
RAVEN (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite that allows for semi-automated reconstruction of genome-scale models. It makes use of published models and/or the KEGG database, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology

If you are using RAVEN in any scientific work, please cite: [R. Agren, et. al, “The RAVEN Toolbox and Its Use for Generating a Genome-scale Metabolic Model for Penicillium chrysogenum,” PLoS Comput. Biol., vol. 9, no. 3, p. e1002980, Mar. 2013.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980).

Please report any technical issues and bugs [here](https://github.com/SysBioChalmers/RAVEN/issues). For other issues, please contact [Eduard Kerkhoven](https://github.com/edkerk).

## Releases
RAVEN can be installed via cloning the Github repository as per below or by downloading and extracting one of the a zipped [release](https://github.com/SysBioChalmers/RAVEN/releases). Please note that the releases do not always represent the most up to date version.

## Installation
If using git please execute the below command in an appropriate directory.

```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```

Don't forget to add RAVEN (the directory extracted from the a release zip file OR the cloned git repository) to you matlab path. In matlab:

```matlab
addpath(genpath('path/to/raven'))
```

### Dependencies:
Required dependencies are currently [libSBML and SBMLToolbox](http://sbml.org/Downloads) used for importing and exporting GEM models in SBML format. Install them and add to the matlab path:

```matlab
addpath(genpath('path/to/libSBML'))
addpath(genpath('path/to/SBMLToolbox'))
```

An optional but recommended dependency is [Gurobi](http://www.gurobi.com/downloads/gurobi-optimizer). Academic licenses are available [here](https://user.gurobi.com/download/licenses/free-academic).

To install the gurobi matlab interface go to the gurobi installation directory and execute `gurobi_setup` as explained [here](https://www.gurobi.com/documentation/6.5/refman/matlab_setting_up_the_guro.html). On Mac OS X in Matlab:

```matlab
cd '/Library/gurobi651/mac64/matlab/'
gurobi_setup
```

To have matlab remember the above installations please save the current path to ’pathdef.m’ *in your MATLAB startup directory*. For instance on Mac OS X:

```matlab
savepath '~/Documents/MATLAB/pathdef.m'
```

*It is highly recommended to set set startup directory of matlab to a directory where you have write access (Preferences->General->start up directory).*

### Setting the solver
Please note that you must now select the solver that you would like to use after you have loaded RAVEN. This can be done like so:

```matlab
setRavenSolver('gurobi')
```

Gurobi is currently the recommended solver (the other available choices are 'mosek' and in the future 'glpk').

## Tutorials
Some tutorials highlighting basic RAVEN functionality can be found in the 'tutorial' folder in the installation directory.

## Pre-trained Hidden Markov Models (HMMs) for KEGG Orthology (KO) protein sets
_For newer RAVEN versions, including GitHub commits after RAVEN 1.9.0_

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