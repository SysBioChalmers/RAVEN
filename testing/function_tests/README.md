# RAVEN function tests

A `matlab.unittest` suite covering the RAVEN source functions, organised as one
test class per source folder (`tQueries`, `tManipulation`, `tSolver`, …).

## Running

```matlab
results = runtests('testing/function_tests');          % whole suite
results = runtests('testing/function_tests/tIO.m');    % one folder
```

A working solver is required for the simulation-based tests; GLPK ships with
RAVEN and is used by default. Some tests additionally need a MILP-capable
solver (Gurobi or SCIP).

## Structure

- `RavenTestCase.m` — abstract base class. Provides the RAVEN root, a fresh
  copy of a small E. coli model before every test, a deterministic solver, and
  helpers (`assumeMILPSolver`, `assumeDependency`, task-model fixtures, figure
  suppression).
- `t<Folder>.m` — one class per source folder; each test method exercises a
  function's main path, favouring structural invariants (counts, membership,
  round-trips, idempotence) over brittle golden values.

## Skipped tests

Tests that need a dependency which may be absent are written with test
*assumptions*, so they are reported as filtered (not failed) when the
dependency is missing. This keeps the suite green in any environment. Skipped
categories: external binaries (BLAST+, DIAMOND, HMMER), a local KEGG dump,
optional solvers (Gurobi/SCIP for MILP), the Optimization Toolbox (`quadprog`),
a Python/pyyaml bridge, and input files for table-driven curation. The full
ftINIT MILP pipeline (in `tINIT`) is guarded on a MILP solver and skipped when
none is available.
