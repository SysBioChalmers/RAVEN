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
| Reaction scoring | `scoreModel` | `scoreComplexModel`, `groupRxnScores` |
| Core MILP | `runINIT` | `ftINITInternalAlg`, scheduled by `getINITSteps` |
| Task gap-filling | `fitTasks` | `ftINITFillGapsForAllTasks` |
| Gene pruning | inline in `getINITModel` | `removeLowScoreGenes` |

What they genuinely share is RAVEN's general machinery: `checkTasks` and
`getEssentialRxns` for task feasibility, `parseTaskList`, `simplifyModel`, the
solver layer, and the model-manipulation and I/O functions.

`fitTasks` has callers outside tINIT, so it stays in `gapfilling/`. `scoreModel`
has none besides `getINITModel`, and is now a wrapper over `scoreComplexModel`
(`omics/`) holding the tINIT argument order, the per-reaction `dataPrecedence`
this method scores with, and its `-Inf` convention for a gene with no data.
