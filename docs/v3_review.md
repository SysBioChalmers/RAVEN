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
  stale branch exists); **FSEOF** is explicitly flagged for a refactor with nicer output;
  new file formats are *not* a priority.

### Already done in `develop3` (not re-proposed here)

So the review doesn't repeat completed work: v3 already removed MetaCyc reconstruction, a
visualization stub, `importExcelModel`, legacy SBML import, Apache POI Excel I/O (replaced
by a dependency-free OOXML writer), `solveQP`, and `combnk`; dropped the Statistics,
Bioinformatics and Text Analytics toolbox dependencies; moved KEGG HMM libraries and the
external binaries out of git to on-demand fetch from the `raven-data` release; reorganized
the repo into functional folders; reformatted docstrings to NumPy style; and added
`progressReport`, `parseRAVENargs` (positional-or-named arguments) and a `matlab.unittest`
suite.

---

## Prioritized to-do list

Tiers are a **recommended** sequencing, not a fixed commitment — see the open questions at
the end. Effort: S ≈ hours–day, M ≈ days, L ≈ week+. Impact is for the RAVEN
reconstruction / strain-design audience.

### Tier 1 — headline v3 work (explicit asks + the biggest leanness wins)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 1 | **MetaNetX / MNXref ID mapping** — revive & refactor the stale `feat/add_MetaNetX` branch into v3 (§1.1) | (b) | M–L | high |
| 2 | **FSEOF redesign** — regression-based target selection, table output, full trajectories (§2) | (c) | M | med–high |
| 3 | **Remove COBRA-duplicate analysis functions** (`runDynamicFBA`, `runPhenotypePhasePlane`, `runProductionEnvelope`, `runRobustnessAnalysis`, `findGeneDeletions`+`qMOMA`) (§4.1) | (a) | S–M | med |
| 4 | **Remove dead/orphaned helpers** (`exportModelToSIF`, `getFullPath`, `getMD5Hash`, `findDuplicateRxns`, `getMetsInComp`, `parseRxnEqu`, `deleteUnusedGenes`) (§4.2–4.3) | (a) | S | med |
| 5 | **Harden RAVEN↔COBRA interop** — round-trip field-loss validator; split `ravenCobraWrapper` (§3) | (b)(d) | M | high |
| 6 | **Central model-field registry** — one declarative table driving add/remove/permute/sort (§5.1) | (d) | L | high |
| 7 | **Verified bug fixes** — `permuteModel` comps, `randomSampling` `showProgress`, others (§7) | — | S | med (correctness) |

### Tier 2 — high-value new capabilities (ambitious; pick a few)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 8 | **Heterologous pathway prediction** over a universal reaction DB (synergizes with #1) (§1.2) | (b) | L | high |
| 9 | **Strain-design upgrade** — EA multi-objective optimization to replace `runSimpleOptKnock` (§1.3) | (b) | M–L | high |
| 10 | **Annotation-coverage / QC report** (MEMOTE-style) (§6.1) | (b) | M | med–high |
| 11 | **Gapfilling family consolidation** (`makeSomething`/`consumeSomething`/`can*`/`checkProduction`) (§5.2) | (d) | M | med |
| 12 | **Deprecate the original INIT path** (`getINITModel`/`runINIT`) in favour of `ftINIT` (§4.4) | (a) | S | med |
| 13 | **Fix HPA omics parsers** broken on current Human Protein Atlas dumps (§6.2) | (c) | M | med |

### Tier 3 — new domains and cheap extras (optional, opportunistic)

| # | Item | Lens | Effort | Impact |
|---|------|------|--------|--------|
| 14 | **OptCouple** — joint knockout + insertion + medium growth-coupling design (§1.3) | (b) | M | med–high |
| 15 | **Regulatory integration** — PROM / CoRegFlux (new application domain) (§1.4) | (b) | M each | med |
| 16 | **Community modeling** — SteadyCom (new application domain) (§1.5) | (b) | M | med |
| 17 | **ROOM / lMOMA** next to `qMOMA`; **E-Flux** next to ftINIT (§1.6) | (b)(c) | S each | low–med |
| 18 | **Smaller consolidations** — deltaG CSV pair, `parseHPA` duplication, WoLF/`parseScores`, KEGG cache plumbing (§5.3–5.5) | (d) | S–L | low–med |

---

## 1. New capabilities

### 1.1 MetaNetX / MNXref ID mapping  *(explicit ask — revive the stale branch)*

**The gap.** Nothing in `annotation/` resolves model metabolite/reaction IDs across
namespaces (BiGG ↔ KEGG ↔ ChEBI ↔ MetaCyc ↔ Rhea ↔ SEED ↔ MNX). `editMiriam`/`extractMiriam`
can *store* any namespace but cannot *map* between them. This is the natural hub for
cross-database interoperability (and a prerequisite for cross-namespace model comparison and
for pathway prediction, #8).

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
  objective; pairs naturally with #8 (insertion candidates come from the pathway DB).
  **Effort M, impact med–high.**

### 1.4 Regulatory–metabolic integration *(new application domain — optional)*

MEWpy's GERM family couples a transcriptional regulatory network to the GEM: **PROM**
(probabilistic) and **CoRegFlux** (regression-based) are the highest-value, lowest-effort
entry points; full Boolean RFBA/SRFBA is more work. RAVEN has no regulatory layer at all, so
this opens an entirely new class of predictions (TF/regulator perturbation phenotypes).
**Effort M per method, impact med.** Pick this if the target users lean metabolic-engineering.

### 1.5 Community / microbiome modeling *(new application domain — optional)*

**SteadyCom** (ReFramed) predicts steady-state species abundances in a stable community; it's
an iterative-LP reformulation that scales well. Optionally add SMETANA-style interaction
scores (MIP/MRO). Community GEMs are a fast-growing area RAVEN can't touch today. **Effort M,
impact med.** Pick this if the target users lean microbiome.

### 1.6 Cheap "while you're in there" additions

- **ROOM / lMOMA** next to the existing `qMOMA` — qualitatively different (often better for
  knockouts) flux predictions. lMOMA is an LP, ROOM a compact MILP. *In COBRA already*, so
  standalone value is modest, but small. **Effort S each.**
- **E-Flux** — continuous expression-scaled bounds for per-condition flux prediction on the
  *full* model, complementing ftINIT's *subnetwork-extraction* approach. **Effort S.**

> **Deliberately excluded:** CAFBA, TFA/TVA (needs a large ΔG dataset), GIMME/iMAT (in
> COBRApy; conceptually covered by RAVEN's INIT family), DifferentialFVA (thin layer over
> FVA). **GECKO-style enzyme constraints are excluded as a *new* capability** because GECKO is
> a sister MATLAB tool in the same ecosystem — the right move is integration, not
> reimplementation.

---

## 2. FSEOF redesign  *(explicit ask)*

`analysis/FSEOF.m` works but is crude. Concrete problems:

1. **Wrong baseline.** Direction/target classification keys off `fseof.results(:,1)`, but
   iteration 1 already enforces `lb = targetMax/iterations` — an already-perturbed state, not
   zero enforcement. There is no true unconstrained baseline.
2. **All-or-nothing, resurrectable monotonicity.** A target is dropped the moment one step is
   non-monotonic, but the flag can be revived in a later iteration — so the computed predicate
   isn't the intended "monotone across the whole scan".
3. **2-point slope.** The slope is a secant between the first and last columns, ignoring all
   intermediate points; noisy and divided by an awkward denominator.
4. **Hand-rolled output.** A 7-column tab-delimited table built with `fopen`/`fprintf`,
   columns mislabeled "Enzyme ID/Name" for reactions, the slope string recomputed twice, and a
   `targets` struct that carries only `.logical`/`.slope` (no trajectories, no direction, no
   fit quality).

**Redesign sketch:**
- Compute an explicit **baseline** at `lb = 0`, then scan `k = 1..iterations` enforcing
  `lb = k/iterations · targetMax`. Store the full `rxns × (iterations+1)` flux matrix `F`
  (baseline in column 1).
- **Target selection via regression**, not step flags: fit each reaction's flux vs. enforced
  target flux (`polyfit` degree 1) over the whole scan. A reaction is a target if its slope is
  significantly non-zero and monotone in one sign (`all(diff(F(j,:)) >= -tol)` or `<= tol`) — a
  single clear predicate. Report **slope = the regression coefficient** plus **R²** so targets
  can be ranked by fit quality.
- Classify and report **amplification, knock-down and knock-out** targets as first-class
  outputs (knock-down/knock-out — reactions driven toward zero — are arguably the most
  actionable and are silently dropped today).
- **Output a MATLAB `table`** (`targetRxn`, `rxnName`, `slope`, `R2`, `direction`,
  `subSystem`, `grRule`), sortable/filterable, with optional `writetable` for the text file.
  Return the full `F` trajectory in the struct for plotting/post-processing.

This overlaps with the already-catalogued back-port notes **FS1/FS2/FS4** (see §8) — fold
them in.

---

## 3. COBRA Toolbox interoperability

`conversion/ravenCobraWrapper.m` converts RAVEN↔COBRA but:

- **No round-trip validator.** The help text *lists* the fields lost in a `raven→cobra→raven`
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

> v3 already breaks backward compatibility, so breaking removals are acceptable. The one
> caveat: several §4.1 functions operate directly on RAVEN structs; removing them tells users
> to convert to COBRA for that analysis, which makes the §3 round-trip validator a soft
> prerequisite. **Confirm appetite** (open question Q2).

### 4.1 COBRA-duplicate analysis functions

All five have **no internal callers** (only tutorials/tests) and are labeled or structured as
COBRA ports:

- **`analysis/runDynamicFBA.m`** — direct port of COBRA `dynamicFBA` (static-optimization
  dFBA); header says "Modified from COBRA Toolbox". No RAVEN-specific value, hard-coded plot.
- **`analysis/runPhenotypePhasePlane.m`** — port of COBRA `phenotypePhasePlane`; calls
  `close all force` (destroys the user's unrelated figures), opens 3 hard-coded figures,
  swallows all errors in a bare `try/end`.
- **`analysis/runProductionEnvelope.m`** + **`runRobustnessAnalysis.m`** — both "Modified from
  COBRA Toolbox"; thin single/double-axis scans over `solveLP`.
- **`analysis/findGeneDeletions.m`** + **`solver/qMOMA.m`** — duplicate COBRA single/double
  gene deletion; the over-expression modes depend entirely on `qMOMA`, which pulls in the
  **Optimization Toolbox** (`quadprog`) the v3 solver layer otherwise avoids and whose own
  comment admits it "never converges good enough" (it hard-codes success). Removing both drops
  a toolbox dependency.

**Proposal:** remove all five (+ `qMOMA`), with their test/tutorial references. If any one is
considered teaching-valuable, demote it and strip the destructive plotting. **Effort S–M,
impact med** (meaningful leanness; aligns with "don't replicate COBRA").

### 4.2 Dead / orphaned I/O helpers

Each has callers only in its own test:

- **`io/exportModelToSIF.m`** — Cytoscape SIF export; niche legacy network format.
- **`io/getFullPath.m`** — 345-line third-party Windows long-path canonicalizer (MATLAB-6.5-era
  fallbacks). Base MATLAB covers any real need.
- **`io/getMD5Hash.m`** — per-OS shell-out (`certutil`/`md5sum`/`md5`); the only production use
  is `getBlast.m` caching DB hashes, and it has a bug (a stray trailing `"` in the `certutil`
  command). Replace that one call site with base-MATLAB `Simulink.getFileChecksum(f,'MD5')` and
  delete the wrapper.

**Verified still needed (keep):** `cleanSheet` (used by `parseTaskList`/`writeExcel`/
`exportToExcelFormat`), `checkFileExistence` (used widely), `makeFakeBlastStructure`.

### 4.3 Redundant query / manipulation helpers

- **`manipulation/findDuplicateRxns.m`** — O(n²) reimplementation of duplicate detection that
  `contractModel` / `simplifyModel(deleteDuplicates)` already do via `unique(S','rows')`. No
  internal callers. Remove; point users to `contractModel`.
- **`queries/getMetsInComp.m`** — thin one-liner used only by tests; has a latent
  error-on-error in its error path. Fold into `getIndexes` (add a `'comps'` type) or remove.
- **`queries/parseRxnEqu.m`** — duplicates parsing already inside `constructS` (its only
  caller), which already returns `mets`.
- **`manipulation/deleteUnusedGenes.m`** — byte-for-byte the same field-trimming as
  `removeReactions(model, [], false, true)` plus a print. Make it a 3-line wrapper or remove.

### 4.4 Deprecate the original INIT path

`getINITModel` → `runINIT` → `fitTasks` is the original tINIT implementation. **`ftINIT` is a
fully separate, faster reimplementation that does not call `runINIT` at all.** The old path is
referenced only by tests/docs — no production caller. Maintaining both means two MILP
formulations, two gap-fill stacks, two scoring entry points. **Do not hard-remove now**
(published Human-GEM-derived models were built with it and users reproduce them): mark
`getINITModel`/`runINIT` **deprecated in v3** with a clear pointer to `ftINIT`, and schedule
removal for v3.x. **Effort S to deprecate.**

---

## 5. Simplifications & internal refactors

### 5.1 Central model-field registry  *(biggest leanness win)*

The dominant smell across `manipulation/` is **hand-maintained model-field ladders**: ~6
functions (`removeReactions`, `permuteModel`, `changeRxns`, `addRxns`, `addMets`,
`addExchangeRxns`, `addTransport`) each carry a near-identical but subtly **divergent** list of
the ~25 rxn / ~14 met / gene / comp fields, copied and drifting apart (e.g. `removeReactions`
handles `spontaneous`/`rxnScores`; `permuteModel` doesn't; `changeRxns` copies `pwys` that
`addRxns` doesn't recognize). The drift has already produced at least one silent correctness
bug (§7, `permuteModel` comps).

**Proposal:** a single declarative field-registry (each field tagged rxn-/met-/gene-/comp-
indexed) consumed by all add/remove/permute/sort/convert operations. This removes the most
code in the repo, fixes the drift bugs structurally, and makes `changeRxns` /
`addRxnsGenesMets` thin wrappers over a shared core. **Effort L, impact high.** The single
highest-leverage refactor for v3.

### 5.2 Gapfilling family consolidation

`gapfilling/` has heavy overlap: `canProduce`/`canConsume` are 2-line wrappers
(`addExchangeRxns` + `haveFlux`); `checkProduction` and `checkRxn` re-implement the same
"add exchange + solve"; **`makeSomething` and `consumeSomething` are ~190-line near-twins**
differing only in the sign of a fake exchange. **Proposal:** one internal
`canExchangeMets(model, direction, mets, …)` behind `canProduce`/`canConsume`/`checkRxn`, and
merge `makeSomething`/`consumeSomething` into one `findGratuitousFlux(model, direction, …)`.
Removes ~250 duplicated lines. The duplicated annotation-field padding both twins use to keep
`solveLP`/`getMinNrFluxes` happy should move into the solver layer (make `solveLP` tolerant of
missing annotation fields) rather than living in callers. **Effort M, impact med.**

### 5.3 Split `ravenCobraWrapper`

See §3 — split the monolith into directional internals.

### 5.4 KEGG `get*FromKEGG` cache plumbing

`getMetsFromKEGG` (288 lines), `getRxnsFromKEGG` (516), `getGenesFromKEGG` (409) triplicate the
"load `kegg*.mat` cache, else parse raw flat files" plumbing. Factor the load-or-build wrapper
into one helper parameterized by entity type; leave the entity-specific parsers. **Effort L,
impact med — lower priority** (core reconstruction, higher risk).

### 5.5 Other small consolidations

- **`annotation/loadDeltaGfromCSV.m` + `saveDeltaGtoCSV.m`** — mirrored load/save pair, four
  near-identical met/rxn blocks; collapse to one `deltaGCSV(model, direction, …)`. **S.**
- **`omics/parseHPA.m` vs `parseHPArna.m`** — extract the shared file-read + header-validation
  into a private helper (combine with the §6.2 fix — same code region). **M.**
- **`localization/getWoLFScores.m` vs `parseScores.m`** — score normalization is reimplemented
  in both; have `getWoLFScores` delegate parsing to `parseScores(file,'wolf')`. Also reconsider
  whether the Linux+Perl-only WoLF runner still earns its place. **S–M.**

---

## 6. Quality / QC additions

### 6.1 Annotation-coverage / QC report

`annotation/` can edit MIRIAMs but cannot *report* annotation completeness. Add
`getAnnotationCoverage` summarizing per-field namespace coverage (what fraction of
mets/rxns/genes carry which namespaces) and flagging unannotated/duplicate-ID entities —
MEMOTE-style, cheap given `extractMiriam` already unpacks the data. Pairs with the MNX mapper
(#1) as the cross-DB annotation hub. **Effort M, impact med–high.**

Relatedly, the catalogued **V2** back-port (see §8) — have `checkModelStruct` return a struct
array of issues (`category`/`target`/`message`) instead of only printing/throwing — would make
model QC programmatically consumable.

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
  `Dbxref`, so the table can feed the MNX/MIRIAM annotation path. **S–M.**

---

## 7. Bug fixes found during the review (do regardless of tiering)

- **`manipulation/permuteModel.m:166–172` (verified).** The `'comps'` case writes remapped
  `rxnComps`/`geneComps` to `model.*` (the input) instead of `newModel.*` (the output); since
  `newModel` was copied earlier, the remap is **silently discarded** — only `metComps`
  survives. One-line fix per field. *Correctness.*
- **`analysis/randomSampling.m:110` (verified).** References `showProgress`, which is never
  defined; errors whenever `sol.f == 0`. Replace with a real flag or drop the branch.
- **`io/getMD5Hash.m` (if kept).** Stray trailing `"` in the `certutil` command string.
- **`analysis/runPhenotypePhasePlane.m:45` (if kept).** `close all force` destroys the user's
  unrelated figures — remove.
- **`optimizeProb.m`** has a dead `if isfield(res,{'dual','rcost'})` block that assigns fields
  to themselves — no-op, remove.
- **Local stray:** `solver/optimizeProb.asv` exists in the working tree but is **not**
  git-tracked — a local autosave artifact; add `*.asv` to `.gitignore` and delete locally.

---

## 8. Relationship to the existing back-port notes

Many granular bug/improvement items are **already catalogued** in the raven-toolbox planning
docs (`docs/reference/matlab_raven_backports.md`, `improvements.md`), targeted at MATLAB RAVEN
from the Python port. Several findings above overlap and should simply be folded into the v3
work rather than re-derived:

- **FSEOF** §2 ⇄ **FS1/FS2/FS4** (regression slope, knockdown/knockout outputs).
- **`randomSampling`** ⇄ **SAMP1/SAMP2/SAMP4** (FVA-based `goodRxns`, `seed`/`nObjectives`
  params, run goodRxns before the inf-bound replacement — the docstring's "3 objectives"
  claim is a known transcription bug).
- **`reporterMetabolites`** ⇄ **RM1** (exact closed-form background correction replacing the
  100k-sample Monte Carlo — faster *and* deterministic).
- **INIT/ftINIT** ⇄ **I4 / FT3 / FT9 / FT8** (per-reaction big-M instead of hard-coded 1000;
  clamp essential forcing to bounds; tree-based `removeLowScoreGenes`).
- **`getModelFromHomology`** ⇄ **H1/H2/H4**; **KEGG pipeline** ⇄ **K1/K2/K3/K15**;
  **`getIndexes`** ⇄ **G1/G5**; **`getElementalBalance`** ⇄ **B2**; **`checkModelStruct`** ⇄
  **V2**; **`addRxns`** ⇄ **A1** (string keyword instead of `eqnType 1/2/3`).

This review was done independently of the Python port (per the brief), but the overlap is
worth exploiting: those notes already contain concrete MATLAB patch shapes.

---

## Open questions for the maintainer

1. **Which Tier-2/Tier-3 capabilities are in scope for v3?** The "ambitious" appetite supports
   several, but #8 (pathway prediction), #9 (EA strain design) and the new-domain items
   (#15 regulatory, #16 community) are each substantial. A reasonable headline set:
   **#1 MetaNetX + #8 pathway prediction + #9 strain design** (one coherent reconstruction →
   design story), with regulatory/community deferred unless they match the target users.
2. **Removal appetite for the COBRA-duplicate analysis functions (§4.1).** Remove outright, or
   keep one or two for teaching and just fix their hygiene? This decision sets whether the §3
   round-trip validator is a hard prerequisite.
3. **Where should new strain-design / community code that needs metaheuristics live**, given
   the dependency question — lean on the Global Optimization Toolbox (`gamultiobj`), or keep
   v3's "minimize toolbox dependencies" stance and ship a self-contained EA?
4. **MNXref data hosting** — confirm the `raven-data` release is the right home for the MNXref
   4.5 assets (mirroring KEGG HMMs/binaries), and whether to ship a prebuilt `.mat` (7 MB,
   fast) or only the gzipped TSVs (slower first build, always current).
