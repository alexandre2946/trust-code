@page Valid_WriteReport How to write a validation report

This page describes how to author a TRUST (or Baltik) validation report.
It focuses on **test-case declaration** through the `trustutils.run` API.
Postprocessing and plotting (`trustutils.plot`, `trustutils.visit`,
`trustutils.widget`) are covered on a separate page.

See also @ref Valid_RunReport for running an existing report, and
@ref Dev_ValidationFormAPI for the full doxygen API reference of `trustutils`.

[TOC]

## What is a validation report                              {#valid-write-what}

A validation report is a **Jupyter notebook** paired with a **`src/`** directory,
sitting under `Validation/Rapports_automatiques/<Category>/<FormName>/`. It plays
three roles at once:

- a human-readable narrative explaining what is being validated and why,
- a runnable artifact that can be exported to PDF or HTML,
- the source from which **non-regression (NR) tests** are automatically derived.

The intended author is a TRUST or Baltik developer adding or qualifying a
feature. The intended consumer is a future reader who needs to understand what
was tested, plus the validation harness, which re-executes the form regularly.

The sample form `Validation/Rapports_automatiques/Verification/SampleFormJupyter`
is a living example of every feature mentioned on this page.


## Anatomy of a form on disk                                {#valid-write-anatomy}

A form has a small, fixed shape:

```
MyForm/
├── MyForm.ipynb         # the report itself
└── src/
    ├── *.data           # template data files (with $-keys for substitution)
    ├── *.med / *.geo    # meshes and other inputs
    └── python_modules/  # reusable Python (see below)
```

Two invariants:

- the notebook is invoked from the form root (where `src/` lives),
- everything under `src/` — and only that — is copied into `build/` when the
  form runs. `build/` is regenerated on demand and is not under version control.

Files such as legacy `src/prepare`, `src/pre_run`, `src/post_run` shell scripts
and `src/liste_cas_exclu_nr` may still appear in older forms. New forms should
not use them; see [Patterns to avoid](#valid-write-legacy).


## A minimal first form                                     {#valid-write-minimal}

Five cells are enough for a runnable form. Markdown cells are added around them
to describe the study.

```python
# Cell 1 - boilerplate
from trustutils import run

run.introduction("Author Name", "21/11/2025")
run.TRUST_parameters()
```

```python
# Cell 2 - reset state and declare cases
run.reset()
run.addCase(".", "my_case.data")
```

```python
# Cell 3 - run
run.runCases()
```

```python
# Cell 4 - performance summary
run.tablePerf()
```

The form is then opened in JupyterLab and executed top to bottom. PDF/HTML
export and full NR execution are driven by `Run_fiche` (see
@ref Valid_RunReport).

Use `from trustutils import run`. The longer `from trustutils.jupyter import run`
still works in many existing forms but is no longer the recommended path.


## Declaring test cases                                     {#valid-write-decl}

Cases are declared **before** running them, all in one cell (or a few),
through the `run` module. The first call to any `addCase*` function
materialises `build/` by copying `src/` into it.

> **Gotcha.** Do **not** `import matplotlib.*` or `from trustutils import plot`
> in any cell that runs before `run.runCases()`. The non-regression
> extraction process executes every cell up to that point in a headless
> context, and a matplotlib import there breaks extraction. Keep all
> plotting imports in postprocessing cells. See @ref valid-post-mpl-order.


### `run.addCase` — the default                             {#valid-write-addcase}

```python
run.addCase("Couette", "Couette_cylindrique.data", nbProcs=2)
run.addCase(".",       "diffusion.data")
```

Use this for the vast majority of cases. The first argument is the directory
**relative to `src/`** (use `"."` for the root). Common keyword arguments:

- `nbProcs` (default `1`) — number of MPI ranks for the case,
- `execOptions` — extra command-line options forwarded to the `trust` binary,
- `excluNR=True` — exclude this case from non-regression test extraction
  (see [Non-regression integration](#valid-write-nr)),
- `pre_run=`, `post_run=` — Python hook callables (see
  [Hooks](#valid-write-hooks)).

`addCase` returns the underlying `TRUSTCase` object, which can be stored for
later manipulation:

```python
case = run.addCase("Couette", "Couette_cylindrique.data", nbProcs=2)
case.partition()      # required when nbProcs > 1, see below
```


### `run.addCaseFromTemplate` — parameter sweeps             {#valid-write-template}

To generate a family of cases from a single template, write `$keyword`
placeholders in the `.data` file and substitute them at declaration time:

```python
for n in (8, 16, 32):
    run.addCaseFromTemplate(
        "diffusion.data",
        targetDirectory=f"mesh_{n}",
        dic={"number": n},
    )
```

This is the preferred way to express mesh-refinement studies, discretization
sweeps, or scheme comparisons.

Pass `targetData="new_name.data"` if the generated file should have a different
name than the template. All the usual kwargs (`nbProcs`, `excluNR`, `pre_run`,
`post_run`) are accepted.


### `TRUSTCase.copy(...).substitute(...)` — the escape hatch  {#valid-write-escape}

When parameter substitution is not enough — for instance to derive a case from
a previously derived case, or to manipulate the file tree programmatically —
use the lower-level methods on `TRUSTCase`:

```python
refined = run.addCaseFromTemplate("diffusion.data", "fine", {"number": 32})
medium  = refined.copy("NewData.data", targetDirectory="medium")
medium.substitute("Nombre_de_Noeuds 32 32", "Nombre_de_Noeuds 16 16")
run.addCase(medium)
```

Note that `copy()` does **not** carry `nbProcs`, `excluNR`, `pre_run` or
`post_run` over — these must be re-specified on the copy.


### Parallel cases — `nbProcs` and `.partition()`            {#valid-write-parallel}

A case declared with `nbProcs > 1` must be partitioned before run. Partitioning
calls `trust -partition` internally:

```python
case = run.addCase("Couette", "Couette_cylindrique.data", nbProcs=2)
case.partition()                              # cuts the domain in N pieces
# or
case.partition(overwritePartition=False)      # keep the Partition block from the dataset
```

By default (`overwritePartition=True`) the partition declared in the dataset is
overwritten, and the domain is cut along a single direction — the original
historical behaviour. Pass `overwritePartition=False` to reuse the
user-specified `Partition` block from the `.data` file.


## Setup work — do it in the notebook                       {#valid-write-setup}

Any preparation that historically lived in a `src/prepare` bash script
(generating meshes, converting MED files, sourcing environments, deriving
intermediate datasets, ...) should now be a regular Python cell in the
notebook, executed **before** the `addCase*` calls.

```python
run.reset()
run.initBuildDirectory()   # materialise build/ explicitly if needed

# --- preparation work, plain Python ---
run.useMEDCoupling()       # makes 'import medcoupling' work
import medcoupling as mc
# ... generate meshes, write files into BUILD_DIRECTORY, etc.

# --- then declare the cases ---
run.addCaseFromTemplate("diffusion.data", "coarse", {"number": 8})
```

Helpful entry points:

- `run.BUILD_DIRECTORY` — absolute path of the build directory,
- `run.initBuildDirectory()` — explicitly copy `src/` into `build/`
  (otherwise this happens implicitly on the first `addCase*`),
- `run.useMEDCoupling()`, `run.useLataTools()` — load the corresponding
  Python bindings into `sys.path`.

If the preparation logic grows past a few lines, move it into
`src/python_modules/` and import it from the notebook
(see [Reusable Python](#valid-write-pymods)).


## `pre_run` and `post_run` hooks                            {#valid-write-hooks}

Each case can be assigned **two Python callables** that run immediately before
and after the `trust` execution. They are the modern replacement for the
legacy `src/pre_run` and `src/post_run` bash scripts, and they are what makes
non-trivial validation patterns (restarts in particular) integrate cleanly
with the non-regression machinery.


### Why hooks exist                                          {#valid-write-hooks-why}

Some validation scenarios cannot be expressed as a single `trust` invocation.
The canonical example is a **restart**: case `B` needs a `.sauv` file produced
by case `A`. Doing this with shell scripts entangles the form with bash
plumbing that the NR extraction process cannot easily reason about, and that
typically resorts to manual file copies and hard-coded paths.

Python hooks solve this. A `pre_run` is just a function: it can call any
TRUST machinery, produce input files, or even launch other `TRUSTCase`s and
wait for them. And because the runtime knows it is a Python callable, it can
**execute it in NR-extraction mode too**, so that the restart case generated
for non-regression also has its `.sauv` file ready.


### Signature and how to attach them                         {#valid-write-hooks-sig}

```python
def pre_run(case_path, case_name):
    """Called just before 'trust <case_name>.data' is executed."""
    ...

def post_run(case_path, case_name):
    """Called just after the case has finished successfully."""
    ...

run.addCase(".", "my_case.data", pre_run=pre_run, post_run=post_run)
```

- `case_path` is the absolute directory of the case inside `build/`,
- `case_name` is the dataset name without the trailing `.data`.

The same `pre_run=` / `post_run=` keyword arguments are accepted by
`run.addCase`, `run.addCaseFromTemplate`, `TRUSTCase()` and `TRUSTCase.copy()`.
A case cannot have **both** a Python hook and a `src/pre_run` (or
`src/post_run`) bash script of the same kind — the runtime will refuse to run
the case and ask for the bash file to be deleted.


### Lifecycle                                                {#valid-write-hooks-lifecycle}

For each case `run.runCases()` orchestrates:

```
                  +-- pre_run(case_path, case_name)
                  |
   run.runCases() +-- trust <case_name>.data <nbProcs> <execOptions>
                  |
                  +-- post_run(case_path, case_name)
```

In the default **parallel mode**, the `pre_run` runs synchronously, the
`trust` invocation is submitted as a subprocess, and the `post_run` is
attached as a callback that fires once the case finishes. In sequential mode
(`run.runCases(preventConcurrent=True)`) all three steps run strictly in
order on the same Python thread before the next case starts.


### Launching cases from a hook                              {#valid-write-hooks-subcases}

A `pre_run` (or `post_run`) may itself construct a `TRUSTCase` and call
`.run()` on it. The scheduler detects these sub-cases, treats them as
dependencies of the parent case, and waits for them to finish before
launching the parent.

This is the recommended pattern for restart-based validation:

```python
def pre_run(case_path, case_name):
    pre = run.TRUSTCase(".", "diffusion.data").copy(
        case_name.replace("_restart", "_pre_run") + ".data",
        targetDirectory=case_path + "/pre_run",
    )
    pre.substitute_template({"number": 20})
    pre.run()      # blocks (or registers as dependency) until .sauv is produced

run.addCaseFromTemplate(
    "diffusion_restart.data", "k20",
    dic={"number": 20, "sauv": "./pre_run/diffusion_pre_run_20_pb.sauv"},
    targetData="diffusion_restart_20.data",
    pre_run=pre_run,
)
```

The form `Validation/Rapports_automatiques/Verification/TestJupyterValidation/Test_Jupyter_pre_run_with_restart`
is the canonical worked example.


### Behaviour under NR extraction                            {#valid-write-hooks-nr}

When the form is reloaded by the NR-extraction machinery
(`run.isExtractingNR()` is `True`), the runtime calls registered Python
`pre_run`s — so that the artifacts they produce (typically restart files)
exist when the NR test runs — but it does **not** call `post_run`s, and it
does not start the main `trust` invocations from the notebook.

Two practical consequences for hook authors:

- `pre_run` is the right place for anything an NR run also needs,
- avoid putting logic in `post_run` that is required for the case itself to
  be valid; put it in `pre_run` or in the dataset.


### Hooks vs legacy bash scripts                             {#valid-write-hooks-vs-bash}

The historical `src/pre_run` and `src/post_run` shell scripts still work for
back-compatibility, but new forms should always prefer the Python hooks:

- Python hooks participate in NR extraction; bash scripts do not.
- Python hooks participate in dependency tracking when they launch
  sub-cases; bash scripts cannot.
- Python hooks live next to the case declaration in the notebook, which
  keeps the form self-contained and reviewable.

A case may not mix the two: if a `pre_run=` callable is provided, the runtime
expects no `src/pre_run` script (and vice-versa).


## Reusable Python                                          {#valid-write-pymods}

Anything more than a few lines — mesh generation helpers, intricate `pre_run`
bodies, custom postprocessing — belongs in `src/python_modules/foo.py` and is
imported by the notebook:

```python
import my_helpers as mh
mh.generate_mesh(n=32)
```

`src/python_modules/` is automatically on `sys.path` when the form runs.

One constraint: **do not import `trustutils.run` from these modules**. The
orchestration (calls to `reset`, `addCase`, `runCases`, ...) belongs in the
notebook itself, so that the form remains the single place that describes the
test plan. The modules should expose pure helpers operating on files, arrays
or `TRUSTCase` instances passed in by the notebook.


## Non-regression integration                               {#valid-write-nr}

Every case declared in a form automatically becomes an NR test. By default,
NR truncates the run to its first **three timesteps**.

- To change the NR timestep cap for a given case, add the magic comment
  `# NB_PAS_DT_MAX n #` (with `n` an integer) directly in the `.data` file.
- To exclude a case from NR entirely, pass `excluNR=True` when adding it:
  ```python
  run.addCase(".", "long_running_case.data", excluNR=True)
  ```
- Do **not** use `src/liste_cas_exclu_nr`. It is deprecated for Jupyter forms
  and the runtime now raises if it lists tests that are not in the report.

When in doubt about NR behaviour, check `run.isExtractingNR()` from a cell to
branch on it explicitly. Keep the case-declaration cell side-effect-free
beyond what NR tolerates: it will be re-imported in that mode.


## Patterns to avoid                                        {#valid-write-legacy}

Older forms exhibit a number of patterns that we no longer recommend.
Authors should not introduce these in new forms, and should migrate them
opportunistically when touching existing forms.

| Avoid                                              | Use instead                                                              |
|----------------------------------------------------|--------------------------------------------------------------------------|
| `from trustutils.jupyter import run`               | `from trustutils import run`                                             |
| `src/prepare` bash script                          | a plain Python cell in the notebook (or `src/python_modules/`)            |
| `src/pre_run`, `src/post_run` bash scripts         | `pre_run=` / `post_run=` Python callables ([Hooks](#valid-write-hooks))  |
| `src/liste_cas_exclu_nr`                           | `excluNR=True` on the case                                                |
| `run.executeScript("...")` / `run.executeCommand`  | reusable Python in `src/python_modules/`                                  |
| Hand-rolled file-tree copies and `sed` of `.data`  | `run.addCaseFromTemplate(...)` or `TRUSTCase.copy().substitute*()`        |


## Authoring loop                                           {#valid-write-loop}

While authoring, the typical loop is to open the notebook in JupyterLab from
the form root, edit and re-execute cells interactively, and only invoke
`Run_fiche` once for full PDF export or NR extraction (see
@ref Valid_RunReport).

Two switches worth knowing:

- `run.saveFormOutput()` opts in to keeping cell outputs in the saved
  notebook (off by default to keep diffs clean),
- `run.runCases(verbose=True)` prints the `trust` console output to the
  notebook for debugging; `run.runCases(preventConcurrent=True)` forces
  strictly sequential execution.


## Pointers                                                 {#valid-write-pointers}

- @ref Valid_PostprocessReport — postprocessing and presenting results.
- @ref Valid_RunReport — running and exporting an existing report.
- @ref Dev_ValidationFormAPI — full doxygen API of the `trustutils` package.
- `Validation/Rapports_automatiques/Verification/SampleFormJupyter` — showcase
  form covering every feature.
- `Validation/Rapports_automatiques/Verification/TestJupyterValidation/Test_Jupyter_pre_run_with_restart`
  — canonical restart-via-hook example.
