# RAVEN 3 — development review and roadmap

A review of the RAVEN toolbox aimed at the next major release (**v3**, currently on the
`develop3` branch). The brief: make the codebase **leaner** (remove old/unused/low-value
functions) while gaining **a few critical new capabilities** that increase value and
usability — with an appetite for ambitious additions.

This document is a proposal and a prioritized to-do list, not a record of work already done.

---

## Scope and method

- Reviewed against the **current `develop3` tree**, not `develop`. The four review lenses:
  **(a)** unnecessary functions to remove, **(b)** missing functionality worth adding,
  **(c)** existing functionality worth expanding, **(d)** functions worth simplifying.
- Method: subsystem-by-subsystem audit of every `.m` file outside `software/`, plus
  ecosystem research on MetaNetX/MNXref and the MEWpy / ReFramed / CAMEO Python tools to
  find capabilities **not** already in COBRA or RAVEN.
- Guard-rails taken from the brief: **don't reimplement COBRA**; interoperability with the
  **COBRA Toolbox (MATLAB)** is wanted; **MetaNetX mapping** is explicitly wanted (an old
  stale branch exists); new file formats are *not* a priority.

### Already done in `develop3`

So the review does not repeat completed work. v3 already:

- Removed MetaCyc reconstruction, a visualization stub, `importExcelModel`, legacy SBML
  import, Apache POI Excel I/O (replaced by a dependency-free OOXML writer), `solveQP`,
  and `combnk`.
- Dropped the Statistics, Bioinformatics and Text Analytics toolbox dependencies.
- Moved KEGG HMM libraries and external binaries out of git to on-demand fetch from the
  `raven-data` release.
- Reorganized the repo into functional folders; reformatted docstrings to NumPy style.
- Added `progressReport`, `parseRAVENargs` (positional-or-named arguments), and a
  `matlab.unittest` suite.
- Fixed verified bugs: `permuteModel` comps remap, `randomSampling` undefined variable,
  `optimizeProb` dead no-op block.
- Backported from raven-toolbox: FSEOF regression-based target selection with knockdown /
  knockout classification (FS1/FS2/FS4); `randomSampling` FVA-based goodRxns and
  reproducible seed (SAMP1/SAMP2/SAMP4); `reporterMetabolites` exact closed-form background
  correction (RM1); ftINIT per-reaction big-M, essential clamping, tree-based gene removal
  (FT3/FT8/FT9); `runINIT` per-reaction big-M (I4); `getModelFromHomology` bidirectional /
  bestHitsOnly / bitscore default / token-based GPR rewrite (H1/H2/H4); KEGG flat-file
  EQUATION-field parsing, undefined-stoichiometry flagging, structured rxn_flags, K15
  cutoff recalibration (K1/K2/K3/K15); `getIndexes` hash-map and islogical fix (G1/G5);
  `getElementalBalance` empty-reaction fix (B2); `checkModelStruct` struct-array issues
  output (V2); `findPotentialErrors` public return + parse-tree detection (S1/S2);
  `addRxns` string keyword for metsBy (A1).
- Removed dead/orphaned I/O helpers: `exportModelToSIF`, `getFullPath`, `getMD5Hash`
  (replaced by inline Java MD5 in `getBlast`).
- Made `dispEM` calls native: replaced all 317 call sites across 67 files with
  `warning()`/`error()` and added `ravenList` for formatted item lists.
- Central model-field registry `ravenModelFields()` implemented (§5.1);
  `permuteModel` and `removeReactions` refactored to use it; `addExchangeRxns`
  and `addTransport` patched for `pwys`/`spontaneous` drift and wrong defaults.
- Removed `getMetsInComp` (thin; callers can use `model.metComps == J` directly);
  collapsed `loadDeltaGfromCSV`/`saveDeltaGtoCSV` into `deltaGCSV(model,direction,…)`;
  rewrote `findDuplicateRxns` with `unique(…,'rows')` (O(n log n) vs O(n²));
  made `deleteUnusedGenes` a thin wrapper over `removeReactions`.
- `getWoLFScores` already delegates to `parseScores` — no change needed.

---

## Prioritized to-do list

Tiers are a **recommended** sequencing, not a fixed commitment — see the open questions at
the end. Effort: S ≈ hours–day, M ≈ days, L ≈ week+. Impact is for the RAVEN
reconstruction / strain-design audience.

### Tier 1 — headline v3 work (explicit asks + the biggest leanness wins)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 1 | **MetaNetX / MNXref ID mapping** — revive & refactor the stale `feat/add_MetaNetX` branch into v3 (§1.1) | (b) | M–L | high |
| 2 | **Harden RAVEN↔COBRA interop** — round-trip field-loss validator; split `ravenCobraWrapper` (§3) | (b)(d) | M | high |

### Tier 2 — high-value new capabilities (ambitious; pick a few)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 6 | **Heterologous pathway prediction** over a universal reaction DB (synergizes with #1) (§1.2) | (b) | L | high |
| 7 | **Strain-design upgrade** — EA multi-objective optimization to replace `runSimpleOptKnock` (§1.3) | (b) | M–L | high |
| 8 | **Annotation-coverage / QC report** (MEMOTE-style) (§6.1) | (b) | M | med–high |
| 9 | **Gapfilling family consolidation** (`makeSomething`/`consumeSomething`/`can*`/`checkProduction`) (§5.2) | (d) | M | med |
| 10 | **Deprecate the original INIT path** (`getINITModel`/`runINIT`) in favour of `ftINIT` (§4.4) | (a) | S | med |
| 11 | **Fix HPA omics parsers** broken on current Human Protein Atlas dumps (§6.2) | (c) | M | med |

### Tier 3 — new domains and cheap extras (optional, opportunistic)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 12 | **OptCouple** — joint knockout + insertion + medium growth-coupling design (§1.3) | (b) | M | med–high |
| 13 | **Regulatory integration** — PROM / CoRegFlux (new application domain) (§1.4) | (b) | M each | med |
| 14 | **Community modeling** — SteadyCom (new application domain) (§1.5) | (b) | M | med |
| 15 | **ROOM / lMOMA** next to `qMOMA`; **E-Flux** next to ftINIT (§1.6) | (b)(c) | S each | low–med |
| 14 | **HPA omics parser consolidation** — extract shared file-read + header-validation into `omics/private/` helper; combine with §6.2 column-name fix (§5.5) | (d) | M | low–med |

---

## 1. New capabilities

### 1.1 MetaNetX / MNXref ID mapping  *(explicit ask — revive the stale branch)*

**The gap.** Nothing in `annotation/` resolves model metabolite/reaction IDs across
namespaces (BiGG ↔ KEGG ↔ ChEBI ↔ MetaCyc ↔ Rhea ↔ SEED ↔ MNX). `editMiriam`/`extractMiriam`
can *store* any namespace but cannot *map* between them. This is the natural hub for
cross-database interoperability (and a prerequisite for cross-namespace model comparison and
for pathway prediction, #6).

**The asset.** The stale branch `origin/feat/add_MetaNetX` (RAVEN 2.8.3 era) already
implements a full MNXref pipeline in `external/MetaNetX/`: `mapToMNX` (entry point) →
`buildMNXref` (parse the MetaNetX TSVs) → `mapModelMets` (map mets by existing xref IDs and
by name) → `mapRxnsViaMets` (map reactions by metabolite-set matching) →
`filterMetMNXIDsViaRxns` → `checkMNXMetsConsistency` (collapse to one MNXID by formula+charge)
→ `addMNXannot` (write `metMetaNetXID`/`rxnMetaNetXID` and back-fill external annotations).
Plus a standalone `mapIDsViaMNXref(type, ids, fromDB, toDB)` ID translator.

**What's stale / blocks a clean merge.**
- **Committed data** — `MNXref.mat` (7 MB) and `chebi.dat` (707 KB) are committed to git, and
  `mapIDsViaMNXref` does `load('MNXref.mat')` off the path. This directly contradicts the v3
  direction (KEGG HMMs and binaries were just moved *out* of git to on-demand fetch).
- **Pinned old version** — defaults to MNXref 3.2 (2020); current is **4.5 (2025-08-13, CC BY
  4.0)**. The TSV layout changed across 4.0/4.2 (added `*_depr.tsv`, `chem_isom.tsv`), and the
  parser does brittle positional column cleanup that will silently break on 4.5.
- **No structural matching** — maps only by existing xref IDs and fuzzy names; ignores the
  **InChIKey** column, which is the modern primary metabolite key.
- **Pre-reorg layout** — lives in `external/MetaNetX/`, carries duplicate `struct_conversion/*`
  that has since diverged on develop3.

**Refactor plan** (one focused PR; core is M, the polish items are S each):
1. Relocate to the develop3 layout (e.g. `reconstruction/metanetx/` or `annotation/metanetx/`),
   drop the duplicated `struct_conversion/*`, rebase on develop3's `convertMiriams`.
2. **Delete the committed data**; host MNXref 4.5 (gzipped TSVs, and ideally a prebuilt `.mat`)
   as a `raven-data` release asset and add a `downloadMNXref.m` modeled exactly on the
   develop3 KEGG-HMM / `downloadRavenBinaries` on-demand pattern (sentinel-file check, cache
   under a gitignored dir).
3. Update parsing to the documented **4.5** columns (header-driven `readtable`, handle the
   `#` comment block and deprecations).
4. Add an **InChIKey** match channel to `mapModelMets` (single biggest accuracy win; aligns
   with current tooling like cobrababel / metanetmap).
5. Clean the public API to two entry points: `mapToMNX(model, …)` (annotate a model) and
   `mapIDsViaMNXref(type, ids, fromDB, toDB)` (pure ID translation, reference passed in or via
   the cached loader — no bare `load`). Normalize emitted prefixes to identifiers.org.
6. Add a unit test that maps a small bundled model and asserts non-zero met/rxn MNX coverage.

**Reference files** (MNXref 4.5, tab-delimited, `https://www.metanetx.org/cgi-bin/mnxget/mnxref/<FILE>`):
`chem_prop.tsv` (incl. InChI/InChIKey/SMILES), `chem_xref.tsv` (largest — the ID map),
`reac_prop.tsv` (compartmentalized equations + EC + balance), `reac_xref.tsv`,
`comp_prop/comp_xref.tsv`, and the `*_depr.tsv` deprecation remaps.

### 1.2 Heterologous pathway prediction *(new — strongest distinctive capability)*

Given a target chemical, enumerate shortest heterologous pathways connecting it to host
metabolism by searching a **universal reaction database**, returning ranked sets of reactions
to add (CAMEO's `PathwayPredictor`). **Not in COBRA, not in RAVEN.** It is the single most
distinctive cell-factory-design feature across MEWpy/ReFramed/CAMEO and fits RAVEN's
reconstruction strengths squarely — and the universal DB is exactly the MNXref reaction set
from #1 (or RAVEN's existing KEGG reaction data). Needs a MILP shortest-pathway formulation
over the universal stoichiometry plus namespace harmonization. **Effort L, impact high.**
Sequence it *after* #1 so it reuses the MNXref plumbing.

### 1.3 Strain-design upgrade *(new — biggest upgrade to existing tooling)*

RAVEN's strain design today is `FSEOF` + `runSimpleOptKnock` (exhaustive small-scale
knockouts). The ecosystem has moved well past this:

- **EA / metaheuristic multi-objective optimization** (MEWpy; CAMEO `OptGene`) — evolutionary
  strain design that co-optimizes knockouts **and** up/down-regulation across reaction / gene /
  enzyme layers with Pareto fronts (BPCY vs yield, etc.). MILP OptKnock scales poorly and
  can't express regulation or non-linear objectives. MATLAB already has `ga`/`gamultiobj`
  (Global Optimization Toolbox), so the engine exists; the work is the problem abstractions
  (reaction/gene knockout & up/down-regulation) and evaluation functions (BPCY, WYIELD).
  **Effort M–L, impact high.** Replaces/extends `runSimpleOptKnock`.
- **OptCouple** (CAMEO) — a single MILP co-designing knockouts + heterologous insertions +
  medium changes for strong growth-coupling. Growth-coupling is the central modern design
  objective; pairs naturally with #6 (insertion candidates come from the pathway DB).
  **Effort M, impact med–high.**

### 1.4 Regulatory-metabolic integration *(new application domain — optional)*

MEWpy's GERM family couples a transcriptional regulatory network to the GEM: **PROM**
(probabilistic) and **CoRegFlux** (regression-based) are the highest-value, lowest-effort
entry points; full Boolean RFBA/SRFBA is more work. RAVEN has no regulatory layer at all, so
this opens an entirely new class of predictions (TF/regulator perturbation phenotypes).
**Effort M per method, impact med.** Pick this if the target users lean metabolic-engineering.

### 1.5 Community / microbiome modeling *(new application domain — optional)*

**SteadyCom** (ReFramed) predicts steady-state species abundances in a stable community; it's
an iterative-LP reformulation that scales well. Optionally add SMETANA-style interaction
scores (MIP/MRO). Community GEMs are a fast-growing area RAVEN cannot touch today. **Effort M,
impact med.** Pick this if the target users lean microbiome.

### 1.6 Cheap "while you're in there" additions

- **ROOM / lMOMA** next to the existing `qMOMA` — qualitatively different (often better for
  knockouts) flux predictions. lMOMA is an LP, ROOM a compact MILP. *In COBRA already*, so
  standalone value is modest, but small. **Effort S each.**
- **E-Flux** — continuous expression-scaled bounds for per-condition flux prediction on the
  *full* model, complementing ftINIT's *subnetwork-extraction* approach. **Effort S.**

> **Deliberately excluded:** CAFBA, TFA/TVA (needs a large deltaG dataset), GIMME/iMAT (in
> COBRApy; conceptually covered by RAVEN's INIT family), DifferentialFVA (thin layer over
> FVA). **GECKO-style enzyme constraints are excluded as a *new* capability** because GECKO is
> a sister MATLAB tool in the same ecosystem — the right move is integration, not
> reimplementation.

---

## 3. COBRA Toolbox interoperability

`conversion/ravenCobraWrapper.m` converts RAVEN to COBRA and back, but:

- **No round-trip validator.** The help text *lists* the fields lost in a `raven->cobra->raven`
  (and the reverse) cycle, but nothing *checks* them. A `compareModelStructFields`-style diff
  that reports which fields/values diverge after a round trip would be high-value, is not COBRA
  duplication, and underwrites the Tier-1 removals in §4.1 (if users are told "convert to COBRA
  for that analysis", the conversion must be demonstrably lossless). **Effort M, impact high.**
- **426-line monolith.** Parsing + both full conversions live in one function whose two halves
  share almost nothing. Split into `ravenToCobra` / `cobraToRaven` internals behind a thin
  dispatcher — easier to test field-by-field. **Effort M, impact med.**
- **Cross-namespace comparison.** Once MNX mapping (#1) lands, a thin `comparison/` helper
  could compare models built in *different* namespaces. Today `compareRxnsGenesMetsComps`
  compares raw ID strings, so two models of the same organism in different namespaces appear to
  share nothing. **Effort M, impact med (depends on #1).**

---

## 4. Leaner code — removals

> v3 already breaks backward compatibility, so breaking removals are acceptable.

### 4.1 COBRA-duplicate analysis functions *(kept for teaching)*

All five have **no internal callers** (only tutorials/tests) and are labeled or structured as
COBRA ports. **Decision: keep all five.** These are simple, self-contained implementations
used for teaching; removing them would require users to install COBRA for basic analyses.
Hygiene items to address instead:

- **`analysis/runPhenotypePhasePlane.m`** — calls `close all force` (destroys unrelated
  figures) and swallows all errors in a bare `try/end`. Strip the destructive figure
  management.
- **`analysis/findGeneDeletions.m`** + **`solver/qMOMA.m`** — `qMOMA` pulls in the
  **Optimization Toolbox** (`quadprog`) and its own comment admits it "never converges good
  enough". Consider a warning that `qMOMA` results are unreliable, or replace with an LP
  approximation (`lMOMA`).

### 4.3 Redundant query / manipulation helpers *(done)*

- **`manipulation/findDuplicateRxns.m`** — kept; rewritten with `unique(S','rows')` for
  O(n log n) instead of the original O(n²) nested loop.
- **`queries/getMetsInComp.m`** — removed. The one-liner `model.metComps == J` (where J is
  from `ismember`) is trivial to write inline; the function had a latent error-on-error bug
  and was only called by tests.
- **`queries/parseRxnEqu.m`** — kept. The doc claimed `constructS` was its only caller, but
  `prepINITModel` also calls it. Remove is deferred until the INIT path is cleaned up (§4.4).
- **`manipulation/deleteUnusedGenes.m`** — refactored to a thin wrapper over
  `removeReactions(model, [], 'removeUnusedGenes', true)`; verbose output retained.

### 4.4 Deprecate the original INIT path

`getINITModel` -> `runINIT` -> `fitTasks` is the original tINIT implementation. **`ftINIT` is a
fully separate, faster reimplementation that does not call `runINIT` at all.** The old path is
referenced only by tests/docs — no production caller. Maintaining both means two MILP
formulations, two gap-fill stacks, two scoring entry points. **Do not hard-remove now**
(published Human-GEM-derived models were built with it and users reproduce them): mark
`getINITModel`/`runINIT` **deprecated in v3** with a clear pointer to `ftINIT`, and schedule
removal for v3.x. **Effort S to deprecate.**

---

## 5. Simplifications & internal refactors

### 5.1 Central model-field registry  *(done)*

`utils/ravenModelFields()` returns a 42-element struct array declaring every 1-D RAVEN model
field with its entity type (`'rxn'/'met'/'gene'/'comp'`) and default value. Fields requiring
2-D indexing (`S`, `rxnGeneMat`, `b`) or compartment-index value-remapping are excluded and
handled explicitly by callers.

Consumers refactored: `permuteModel` (replaced ~80 individual `isfield` blocks with a generic
loop + 4-case switch for 2-D fields, fixing silent drift for `pwys`, `spontaneous`, `metNotes`,
`geneFrom`); `removeReactions` (same for ~22 rxn-field and ~5 gene-field blocks).
`addExchangeRxns` and `addTransport` patched for the `pwys`/`spontaneous` gap; `addTransport`
default for `rxnConfidenceScores` corrected (1→NaN) and `rxnDeltaG` (0→NaN).

Remaining consumers to refactor against the registry: `changeRxns`, `addRxns`, `addMets`,
`addRxnsGenesMets`, `ravenCobraWrapper`. The drift bugs in those functions are low-priority
because the registry now exists as the reference — callers that skip it are visibly incomplete.

### 5.2 Gapfilling family consolidation

`gapfilling/` has heavy overlap: `canProduce`/`canConsume` are 2-line wrappers
(`addExchangeRxns` + `haveFlux`); `checkProduction` and `checkRxn` re-implement the same
"add exchange + solve"; **`makeSomething` and `consumeSomething` are ~190-line near-twins**
differing only in the sign of a fake exchange. **Proposal:** one internal
`canExchangeMets(model, direction, mets, ...)` behind `canProduce`/`canConsume`/`checkRxn`, and
merge `makeSomething`/`consumeSomething` into one `findGratuitousFlux(model, direction, ...)`.
Removes ~250 duplicated lines. The duplicated annotation-field padding both twins use to keep
`solveLP`/`getMinNrFluxes` happy should move into the solver layer (make `solveLP` tolerant of
missing annotation fields) rather than living in callers. **Effort M, impact med.**

### 5.3 Split `ravenCobraWrapper`

See §3 — split the monolith into directional internals.

### 5.5 Other small consolidations

- **`annotation/loadDeltaGfromCSV.m` + `saveDeltaGtoCSV.m`** *(done)* — collapsed into
  `deltaGCSV(model, direction, metCsv, rxnCsv)`. The four near-identical met/rxn blocks are
  handled by two shared local functions (`loadField`/`saveField`).
- **`omics/parseHPA.m` vs `parseHPArna.m`** — shared file-open + header-validation loop
  should be extracted into a private helper in `omics/private/`. Address together with the §6.2
  column-name fix (same code region). **M.**
- **`localization/getWoLFScores.m` vs `parseScores.m`** *(already done)* — `getWoLFScores`
  already delegates to `parseScores(file,'wolf')`; no normalization duplication exists.
  Remaining question: whether the Linux+Perl-only WoLF runner still earns its place.

---

## 6. Quality / QC additions

### 6.1 Annotation-coverage / QC report

`annotation/` can edit MIRIAMs but cannot *report* annotation completeness. Add
`getAnnotationCoverage` summarizing per-field namespace coverage (what fraction of
mets/rxns/genes carry which namespaces) and flagging unannotated/duplicate-ID entities —
MEMOTE-style, cheap given `extractMiriam` already unpacks the data. Pairs with the MNX mapper
(#1) as the cross-DB annotation hub. **Effort M, impact med-high.**

### 6.2 Fix HPA omics parsers (currently broken)

`omics/parseHPArna.m` hard-codes HPA v18/v19 column layouts (and `parseHPA.m` v17/v19); current
Human Protein Atlas dumps are well past these, so the functions **silently error on today's
data**. Make column detection header-driven (match by name, not fixed position). **Effort M,
impact med.**

### 6.3 Smaller expansions worth bundling

- **`annotation/editMiriam.m`** — add a `'remove'` mode (today it can only replace/fill/add)
  and optional identifiers.org prefix validation. **S.**
- **`queries/getExchangeRxns.m`** — add an optional output of the exchanged metabolite per
  reaction (callers like `setExchangeBounds` re-derive it from `S`). **S.**
- **`curation/getGeneData.m`** — also surface UniProt/RefSeq xrefs already present in the GFF
  `Dbxref`, so the table can feed the MNX/MIRIAM annotation path. **S-M.**

---

## Open questions for the maintainer

1. **Which Tier-2/Tier-3 capabilities are in scope for v3?** The "ambitious" appetite supports
   several, but #6 (pathway prediction), #7 (EA strain design) and the new-domain items
   (#13 regulatory, #14 community) are each substantial. A reasonable headline set:
   **#1 MetaNetX + #6 pathway prediction + #7 strain design** (one coherent reconstruction ->
   design story), with regulatory/community deferred unless they match the target users.
2. **Where should new strain-design / community code that needs metaheuristics live**, given
   the dependency question — lean on the Global Optimization Toolbox (`gamultiobj`), or keep
   v3's "minimize toolbox dependencies" stance and ship a self-contained EA?
4. **MNXref data hosting** — confirm the `raven-data` release is the right home for the MNXref
   4.5 assets (mirroring KEGG HMMs/binaries), and whether to ship a prebuilt `.mat` (7 MB,
   fast) or only the gzipped TSVs (slower first build, always current).
