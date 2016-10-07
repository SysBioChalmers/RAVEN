# RAVEN
RAVEN (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite that allows for semi-automated reconstruction of genome-scale models. It makes use of published models and/or the KEGG database, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology

If you are using RAVEN in any scientific work, please cite: [R. Agren, et. al, “The RAVEN Toolbox and Its Use for Generating a Genome-scale Metabolic Model for Penicillium chrysogenum,” PLoS Comput. Biol., vol. 9, no. 3, p. e1002980, Mar. 2013.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980).

## Installation
RAVEN can be installed via cloning the Github repository as per below or by downloading and extracting one of the a zipped [release](https://github.com/SysBioChalmers/RAVEN/releases). Please note that the releases do not always represent the most up to date version.

```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```

Don't forget to add raven to you matlab path. In matlab:

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

### Setting the solver
Please note that you must now select the solver that you would like to use after you have loaded RAVEN. This can be done like so:

```matlab
setRavenSolver('gurobi')
```

Gurobi is currently the recommended solver (the other available choices are 'mosek' and in the future 'glpk').

## Tutorials
Some tutorials highlighting basic RAVEN functionality can be found in the 'tutorial' folder in the installation directory.

## Links
For more information on software connected to Genome Scale models please visit the Systems Biology [Github page](https://github.com/SysBioChalmers). For information and publications by the Systems Biology department please visit [SysBio](www.sysbio.se).

## License (MIT)
Copyright (c) 2016 SysBio Chalmers

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.