# libSBML MATLAB binding (vendored)

`TranslateSBML_RAVEN` / `OutputSBML_RAVEN` are libSBML's MATLAB binding, renamed so they
cannot clash with an upstream libSBML installation on the same path. `importModel` and
`exportModel` call the MEX files; `exportModel` additionally calls `getDefaultValues` and
`getStructureFieldnames` directly.

**The other `.m` files here look unreferenced but are not.** The MEX files call them by name
through `mexCallMATLAB`, so no static search of the RAVEN sources will find a caller. Deleting
one breaks SBML I/O only at run time. The names embedded in the binaries are:

| Called from the MEX | Reached from those |
| --- | --- |
| `addLevelVersion`, `CheckAndConvert`, `ConvertFormulaToMathML`, `getStructure`, `isoctave`, `isSBML_Model` | `applyUserValidation`, `getDefaultValues`, `getStructureEnum`, `getStructureFieldnames`, `getValueType` |

To re-check after a libSBML upgrade, scan a binary for NUL-delimited symbol names rather than
grepping the `.m` sources.

Upstream helpers that RAVEN removed, and why: `installSBML` (installs upstream's *unrenamed*
`TranslateSBML`/`OutputSBML`, which would then shadow the RAVEN pair), `isEnabled` /
`isFbcEnabled` / `other.m` (all call the unrenamed entry points, so they error here),
`getSBMLDefaultStruct` (reached from nothing), `Contents.m` and the `other.xml` / `test.xml`
fixtures.

The vendored version is recorded in `VERSION.txt`.
