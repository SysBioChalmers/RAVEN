<img src="./RAVEN.png" width="200px">

[![Current release](https://img.shields.io/github/release/SysBioChalmers/RAVEN/all.svg)](https://GitHub.com/SysBioChalmers/RAVEN/releases/)
[![GitHub Discussions](https://img.shields.io/github/discussions-search?query=repo%3Asysbiochalmers%2raven&label=GitHub%20Discussions)](https://github.com/SysBioChalmers/RAVEN/discussions)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.3689518.svg)](https://doi.org/10.5281/zenodo.3689518)
[![MATLAB File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/112330-raven-toolbox)

## About RAVEN

The **RAVEN** (Reconstruction, Analysis and Visualization of Metabolic Networks) Toolbox is a software suite for MATLAB that allows for semi-automated reconstruction of genome-scale models (GEMs). It makes use of published models and/or the KEGG database, coupled with extensive gap-filling and quality control features. The software suite also contains methods for visualizing simulation results and omics data, as well as a range of methods for performing simulations and analyzing the results. The software is a useful tool for system-wide data analysis in a metabolic context and for streamlined reconstruction of metabolic networks based on protein homology.

The latest stable release is **RAVEN 2.11.3**, on the `main` branch.

> **This `develop3` branch is where RAVEN 3 is under active development**, with the upcoming **3.0.0b1** beta as its first release, targeting a final **3.0.0**. RAVEN 3 is a deliberate breaking release: the codebase has been reorganized into functional folders, optional function arguments can now be passed positionally, by name, or both, external binaries and KEGG reference data are fetched on demand instead of bundled in the repository, several MATLAB toolbox dependencies were dropped, and a number of low-value or duplicate functions were removed alongside new capabilities. Until 3.0.0 is tagged, this branch can still change. See the [migration guide](#documentation) for the full list of changes for users upgrading from RAVEN 2.

For context-specific model extraction, use `ftINIT`; `getINITModel` and `runINIT` are the original tINIT implementation, kept for reproducing models built with it.


## Documentation
The information about downloading, installing and developing RAVEN is included in the [Wiki](https://github.com/SysBioChalmers/RAVEN/wiki). The source code documentation is also available [online](http://sysbiochalmers.github.io/RAVEN/doc/).

Full user documentation, including a migration guide for users upgrading from RAVEN 2 to RAVEN 3, is published at [raven-docs](https://raven-docs.readthedocs.io/en/latest/).


## Cite Us
If you use RAVEN in your scientific work, please cite:
> Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J, Kerkhoven EJ. (2018) RAVEN 2.0: A versatile toolbox for metabolic network reconstruction and a case study on _Streptomyces coelicolor_. PLoS Comput Biol 14(10): e1006541. doi:[10.1371/journal.pcbi.1006541](https://doi.org/10.1371/journal.pcbi.1006541).

All releases are also archived in [Zenodo](https://doi.org/10.5281/zenodo.3689518), so you can cite the specific version of RAVEN used in your study

If you use ftINIT in your scientific work, please cite:
> Gustafsson J, Anton M, Roshanzamir F, Jörnsten R, Kerkhoven EJ, Robinson JL, Nielsen J. (2023) Generation and analysis of context-specific genome-scale metabolic models derived from single-cell RNA-Seq data. Proc Natl Acad Sci 120(6): e2217868120. doi:[10.1073/pnas.2217868120](https://doi.org/10.1073/pnas.2217868120)

For crediting supporting work, please cite doi:[10.1002/msb.145122](http://msb.embopress.org/content/10/3/721) (`tInit`); doi:[10.1371/journal.pcbi.1000859](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000859) (`randomsampling`). For crediting RAVEN 1, cite doi:[10.1371/journal.pcbi.1002980](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980). For more details, see [wiki#cite-us](https://github.com/SysBioChalmers/RAVEN/wiki#cite-us).

## Contact Us
Use [GitHub Discussions](https://github.com/SysBioChalmers/RAVEN/discussions) for support, to ask questions or leave comments.

## Contributing

Contributions are always welcome! Please read the [Contributor guidelines](https://github.com/SysBioChalmers/RAVEN/blob/main/.github/CONTRIBUTING.md) to get started. RAVEN 2.x maintenance happens on the `develop` branch; RAVEN 3 development happens on `develop3`.
