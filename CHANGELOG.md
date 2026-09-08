# Changelog

## 3.0.0b1 — 2026-09-07

First beta toward 3.0. Changes below are relative to RAVEN 2.9.3, the last 2.x release:
a reorganized codebase, fewer toolbox dependencies, hardened I/O and MILP algorithms, plus
capabilities backported from the `raven-toolbox` (Python) port.

### Repository & tooling

* Repo reorganized from a flat layout into functional folders (`analysis/`, `annotation/`,
  `biomass/`, `comparison/`, `conversion/`, `gapfilling/`, `io/`, `manipulation/`, `queries/`,
  `solver/`, `utils/`, etc.).
* Function help blocks reformatted to NumPy-style docstrings; design-rationale, version-history
  and cross-language comments trimmed to current-state prose.
* Allow for positional-or-named function arguments via `parseRAVENargs`. Essential arguments
  are positional, non-compulsory arguments are named.
* Dropped the Statistics, Bioinformatics, Text Analytics and Optimization Toolbox dependencies.
* All 317 `dispEM` call sites across 67 files replaced with native `warning()`/`error()`; new
  `ravenList` for formatted item lists (fixed: was silently dropping its own indent).
* `checkRaven` replaces `checkInstallation`, moved to the repo root, faster checks.
* Universal progress reporting (`progressReport`).
* Poor performing visualization functions removed.

### Reconstruction (KEGG & homology)

* KEGG HMM libraries and external binaries moved out of git to on-demand, SHA256-verified fetch
  from the `raven-data` release repo.
* KEGG artifact set updated 116 → 118.
* Legacy KEGG artifact-builder functions removed, moved to raven-toolbox (Python).
* Recalibration of default HMM cutoffs.
* MetaCyc reconstruction removed, as it persistently yielded high false positive annotations.
* `getModelFromHomology`: recalibration of `minLen` default, token-based GPR rewrite (replaces
  a buggy cleanup regex).
* New `assignCompartments`: functionality-constrained compartment-assignment MILP (WIP).

### Curation

* Genome data utilities for homology-based reconstruction: `getGeneData` (gained a `UniProt`
  output column from the GFF `Dbxref` field; fixed crash when a CDS's `Parent` resolves to no
  gene), `downloadGenomeData`, `processProteinFastaFile`.
* `curateModelFromTables`: fixed missing metabolite-support check.

### I/O

* Legacy pre-L3V1/FBCv2 SBML import removed, only L3V1/FBCv2 remains supported.
* `importExcelModel` removed, to discourage the practice of curating models in Excel.
* New `model.ec` export to Excel (ENZYMES/ENZRXNS sheets).
* Apache POI Excel I/O replaced with a dependency-free OOXML writer.
* Changed YAML file format, to better align with cobrapy and have full compatibility with
  raven-toolbox (Python): fixed quoting; multi-line entries; position-dependent defaults;
  `objective_coefficient` import; empty/missing gene sections; compartment annotations;
  EC-code annotation round-trip.
* Support for `rxnDeltaG`/`metDeltaG` in SBML.

### Manipulation

* `addRxns`: new string keyword for `metsBy`; fixed gene ids shredded by `and`/`or`.
* `addTransport`: fixed row-oriented `metNames`; crash on a compartment with no counterpart
  metabolite; `rxnConfidenceScores`/`rxnDeltaG` defaults.
* `mergeCompartments`: fixed inverted unconstrained warning; stale `compMiriams`.
* `mergeModels`: fixed `metFrom` handling; same-batch id-rename collisions.
* `convertToIrrev`: fixed wrong reverse-copy index.
* `setExchangeBounds`: fixed missing multi-exchange metabolites.
* `closeModel`: fixed `b` growth and `metNotes` padding.
* `replaceMets`: fixed three bugs.
* `deleteUnusedGenes` is now a thin wrapper over `removeReactions`.
* `removeLowScoreGenes`: fixed gene field alignment.
* `permuteModel`: fixed compartment remap.
* `simplifyModel`: fixed dropping task constraints during gap-filling.

### Queries

* `checkModelStruct`: struct-array issues output; fixed two bugs; grRules start/end `and`/`or`
  handling.
* `findPotentialErrors` made public, with parse-tree detection; grRules now parsed into a tree,
  also fixing `removeGenes` and `expandModel`.
* `getExchangeRxns`: fixed index return in both branches, and exchange/transport
  misclassification for boundary-metabolite models.
* `getIndexes`: hash-map lookup, `islogical` fix.
* `getElementalBalance`: empty-reaction fix, new charge-balance check.
* `constructS`/`parseFormulas`: fixed output corruption.
* `findDuplicateRxns` rewritten with `unique(…,'rows')`: O(n log n) instead of O(n²).

### Gap-filling

* New LP/MILP/topological gap-filling algorithms: `gapFillFastCore`, `gapFillSwiftCore`,
 `gapFillMILP` and `gapFillTopological`, in addition to the original `fillGaps` function.
* `makeSomething`/`consumeSomething` merged into `findLeakMetabolite(model, direction, …)`
  (originals kept as thin wrappers).

### Context-specific models (ftINIT / tINIT)

* `getINITModel`/`runINIT` deprecated, in favor of `ftINIT` (a separate, faster
  reimplementation that does not call `runINIT`). Code is still available for now, as
  published Human-GEM-derived models were built with tINIT.
* `ftINIT` gained per-reaction big-M, essential-reaction clamping, tree-based gene removal;
  `runINIT` gained per-reaction big-M.
* ftINIT gained `resolveTies`/`proveAbsGap` for deterministic MILP tie-breaking.
* `ftINIT`: fixed metabolomics `spdiags` crash; MILP handling under the free solver; gap-fill
  failures reported as success; inverted `allowExcretion` sense in `ftINITInternalAlg`.
* `ftINIT` task gap-filling now edits the reference model in place instead of copying it.
* `getMinNrFluxes`: fixed crash on all-non-negative scores.
* `sortReactionOrder`: fixed scoring the wrong columns.
* `checkTasks`: task loop parallelized via `parfor`.

### Analysis & flux sampling

* FSEOF regression-based target selection with knockdown/knockout classification.
* `reporterMetabolites` exact closed-form background correction.
* `getFluxZ`: fixed inverted sign in the zero-variance branch.
* New analysis utilities: `modelSummary`, `walkFluxes`, `traceFluxPath`, `compareFluxes`.
* New `getMinimalMedium`: MILP-based minimal medium finder.
* `randomSampling`: FVA-based `goodRxns` selection, reproducible seed parameter; fixed
  undefined variable; per-reaction LPs parallelized via `parfor`; CHRR and ACHR flux samplers
  added.
* `compareRxnsGenesMetsComps` folded into `compareMultipleModels`.
* New `diffModels`: structured model comparison.

### Biomass

* `fitParameters`: fixed every `xRxn`'s upper bound set from the first.
* `guessComposition`: fixed unknown metabolite's own coefficient.

### Localization

* WoLF PSORT localization loader removed, replaced by DeepLoc2, MULocDeep, COMPARTMENTS,
  UniProt (WIP).

### Solver

* `optimizeProb`: fixed dead no-op block; SCIP translation; gurobi results; MILP guard.
