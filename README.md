# RAVEN
RAVEN software to simulate Genome Scale Metabolic models (GEMs) and includes alot of functionality for reconstruction of GEMs based on different kinds of omics data. 

## Installation
RAVEN can be installed via cloning the Github repository as per below or by downloading and extracting one of the a zipped [releases](https://github.com/SysBioChalmers/RAVEN/releases).

```bash
git clone git@github.com:SysBioChalmers/RAVEN.git
```

Don't forget to add raven to you matlab path. In matlab:

```matlab
addpath(genpath('path/to/raven'))
```

### Dependencies:
Required dependencies are currently [libSBML and SBMLToolbox](http://science.sciencemag.org/content/349/6252/1101) used for importing and exporting GEM models in SBML format. Install them and add to the matlab path:

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

## Tutorials
Some tutorials highlighting basic RAVEN functionality can be found in the 'tutorial' folder in the installation directory.

## Links
For more information on software connected to Genome Scale models please visit the Systems Biology [Github page](https://github.com/SysBioChalmers). For information and publications by the Systems Biology department please visit [SysBio](www.sysbio.se).

## License (MIT)
Copyright (c) 2016 SysBio Chalmers

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.