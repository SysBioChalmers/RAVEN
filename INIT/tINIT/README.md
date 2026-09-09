# tINIT

The original tINIT implementation: `getINITModel` drives the reconstruction and
`runINIT` solves the MILP.

It is supported and stays supported. Use it to reproduce models built with it.
For new work use `ftINIT`, in the parent folder, which is the recommended
method for context-specific model extraction.

The two are separate implementations that share no algorithm code. Where they
appear to overlap they do not:

| | tINIT | ftINIT |
|---|---|---|
| Entry point | `getINITModel` | `prepINITModel`, then `ftINIT` |
| Reaction scoring | `scoreModel` (isozyme/complex scoring fixed to `'max'`, `dataPrecedence` `'reaction'`) | `scoreModel` (`INIT/`), `groupRxnScores` |
| Core MILP | `runINIT` | `ftINITInternalAlg`, scheduled by `getINITSteps` |
| Task gap-filling | `fitTasks` | `fitTasks` (`gapFillMode` `'preMerged'`), `ftINITFillGaps` |
| Gene pruning | inline in `getINITModel` | `removeLowScoreGenes` |

What they genuinely share is RAVEN's general machinery: `checkTasks` and
`getEssentialRxns` for task feasibility, `parseTaskList`, `simplifyModel`, the
solver layer, and the model-manipulation and I/O functions.

Reaction scoring and task gap-filling look forked in the table above only
because the two entry points call their shared functions with different
settings, not because two implementations exist. `scoreModel` (`INIT/`) is
one function for both: `getINITModel` calls it with the fixed argument
combination the original tINIT algorithm needs — a single operator for both
`and`/`or` in a grRule, `dataPrecedence` `'reaction'`, and geneScores rewritten
from `NaN` to `-Inf` for a gene with no data — while `ftINIT` calls it with the
general defaults. Likewise `fitTasks` (`gapfilling/`) is the one task
gap-filling loop for both: `getINITModel` uses its default `gapFillMode`
(`'merge'`, backed by `fillGaps`), `ftINIT` passes `'preMerged'` (backed by
`ftINITFillGaps`, since ftINIT's reference model already contains the
sample's own reactions and needs no per-task merge). That MILP-formulation
split — `fillGaps` merging per task vs. `ftINITFillGaps` expecting a
pre-merged model — is the one place the two genuinely differ, not
`getINITModel` vs. `ftINIT` themselves.
