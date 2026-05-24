@page Valid_RunReport How to execute a validation report

This page describes how to run a TRUST (or Baltik) Jupyter validation report
through the `Run_fiche` script. Authoring a new report is covered in
@ref Valid_WriteReport.

> **Scope.** Only the Jupyter-notebook formalism is documented here. The
> historical PRM format is still accepted by `Run_fiche` but is being phased
> out; if you must touch a PRM report, convert it first
> (see [PRM conversion](#valid-run-prm)).
>
> **Heads up.** The `Run_fiche` script has a handful of quirks documented
> under [Gotchas](#valid-run-gotchas). The script is slated for a rewrite;
> these workarounds are temporary.

[TOC]

## Where to run from                                        {#valid-run-where}

`Run_fiche` is always invoked from the **form root** — the directory that
contains the `<FormName>.ipynb` file and its `src/` subdirectory.

- Running it from `src/` is an error and the script will exit.
- The form root must contain exactly one `.ipynb` file.
- All artifacts (datasets after substitution, `.sauv`, `.son`, `.lata`,
  `rapport.pdf`, ...) end up under `build/` next to the notebook.


## Common recipes                                           {#valid-run-recipes}

### Author interactively                                    {#valid-run-recipe-author}

```
Run_fiche
```

Opens the notebook in JupyterLab in your browser. This is the right mode
while writing the form: re-execute cells one by one, inspect intermediate
values, iterate on plots.


### Build a parallel PDF report (recommended)               {#valid-run-recipe-parallel}

```
Run_fiche -parallel_sjob
```

Executes the notebook end-to-end and writes `build/rapport.pdf`. Cases
run concurrently through the **Sserver** local queue (see
[The Sserver](#valid-run-sserver)). This is the right mode for a full run.

> Note that `-parallel_sjob` **already produces the PDF**. Do not add
> `-export_pdf`; the combination has a quirk
> (see [Gotchas](#valid-run-gotchas)).

The non-Sserver alternative is `Run_fiche -parallel_run [N]`, which runs
cases in parallel inside the notebook's Python process using a simple
internal scheduler. Prefer `-parallel_sjob` for non-trivial validation runs:
the Sserver is more robust under recursive case launches (e.g. cases
spawned from `pre_run` hooks).


### Build a sequential PDF report                           {#valid-run-recipe-seq}

```
Run_fiche -export_pdf
```

Executes the notebook end-to-end, sequentially, and writes
`build/rapport.pdf`. Useful for small reports or for debugging an issue you
suspect comes from concurrency.


### Re-render the PDF without rerunning TRUST               {#valid-run-recipe-notrun}

```
Run_fiche -not_run -export_pdf
```

Re-executes the notebook but skips `trust` invocations, Python `pre_run` /
`post_run` hooks, and bash `prepare` / `pre_run` / `post_run` scripts. The
notebook reuses whatever is already in `build/`. The right mode while
iterating on plotting, postprocessing, or narrative text after a full run
has populated `build/`.


### Get just the LaTeX                                       {#valid-run-recipe-tex}

```
Run_fiche -export_tex
```

Same as `-export_pdf` but stops at `build/rapport.tex`. Useful when you want
to tweak the LaTeX manually before running `xelatex`, or to edit the report
title (see the metadata trick at the bottom of `SampleFormJupyter`).


### Bootstrap a notebook from existing data files           {#valid-run-recipe-create}

```
Run_fiche -create
```

Creates a minimal `<FormName>.ipynb` next to `src/`, pre-populated with one
`run.addCase` per `.data` found in `src/`. Useful when starting a new form
from a set of pre-existing datasets.


### Open the showcase form                                  {#valid-run-recipe-doc}

```
Run_fiche -doc
```

Copies `Validation/Rapports_automatiques/Verification/SampleFormJupyter`
to `$TRUST_TMP/` and opens it. Quickest way to see every `trustutils`
feature in action.


## The Sserver                                              {#valid-run-sserver}

The Sserver is a small local job queue shipped with TRUST
(`$TRUST_ROOT/bin/Sjob/`). It throttles concurrent `trust` runs so that a
validation report does not saturate your workstation, and it handles
recursive case launches (cases spawned from `pre_run` hooks) cleanly.

**Start it.** Either let `Run_fiche -parallel_sjob` start it for you (it
will, if `TRUST_WITHOUT_HOST=1` and no Sserver is already running), or
start it explicitly:

```
Sserver <N>      # expose N cores to the queue
```

Re-running `Sserver <N>` against an existing server just adjusts the proc
count.

**Inspect and manage.**

```
Squeue           # list jobs and their status (R/W/S)
Sdelete <id>     # kill a single job and its group
```

**Stop it.** The Sserver is **not** stopped automatically at the end of
`Run_fiche`. Clean up explicitly when you are done with a validation
session:

```
Sserver -1       # stop the server, leave running jobs alone
Sserver -2       # kill all jobs cleanly, then stop the server
```

**Behaviour change with reports.** Whenever an Sserver is alive on the
machine, `trustutils.run` detects it automatically (regardless of flags)
and routes every case through it in parallel mode. A stale Sserver from a
previous session therefore changes how a subsequent `Run_fiche` run
behaves — kill it between sessions to avoid surprises.

See `$TRUST_ROOT/bin/Sjob/README` for the full reference (Salloc options,
advanced client commands).


## Useful side flags                                        {#valid-run-flags}

These can be combined with the PDF-producing recipes above.

- `-o my_report.pdf` — change the PDF output name (default `rapport.pdf`).
- `-xpdf` — open the PDF after the build with the first available viewer.
- `-no_visit` — skip VisIt figures (much faster on heavy VisIt forms, and
  required when VisIt isn't built).
- `-exec /path/to/trust` — override the `trust` binary (default `$exec`).
- `-compare_pdf reference.pdf` — diff the produced PDF against a reference
  PDF after the build. Requires a PDF-producing mode (`-export_pdf`,
  `-parallel_sjob` or `-parallel_run`). The typical use is **checking for
  regressions after a development**: keep a known-good PDF of the report
  on the side, then run
  `Run_fiche -parallel_sjob -compare_pdf good_rapport.pdf` after your
  changes.


## Archiving — required before requesting a PR              {#valid-run-archive}

```
Run_fiche -archive
```

Packages the build directory into a `preserve.tgz` that contains everything
needed to rebuild the report later with a future version of TRUST (the
`.data` files after substitution, mesh files, probes, the `used_files`
manifest, etc.).

`-archive` is **not** part of the day-to-day authoring loop, but every
validation form must support it. Before opening a PR that adds or modifies
a validation report, run:

```
Run_fiche -export_pdf      # or -parallel_sjob, to populate build/
Run_fiche -archive         # must succeed
```

The second invocation reuses `build/`; if it fails, the form is not
PR-ready (typically: a file used during the run was not declared through
the `trustutils.run` API and is missing from `build/used_files`).


## Current gotchas                                          {#valid-run-gotchas}

Until `Run_fiche` is rewritten, these quirks are worth knowing.

### `-export_pdf` silently disables `-parallel_*`            {#valid-run-gotcha-parallel-pdf}

`Run_fiche -parallel_run -export_pdf` runs **sequentially**. The
`-export_pdf` branch of the script does not forward `-parallel_*` flags to
the notebook's runtime.

`Run_fiche -parallel_sjob -export_pdf` *appears* to work, but only
incidentally: parsing `-parallel_sjob` autostarts the Sserver as a side
effect, and `trustutils.run` then auto-detects the running Sserver and
switches to parallel mode on its own. Fragile.

**Recommendation.** Don't combine `-export_pdf` with anything `-parallel_*`.
For a parallel PDF report, use `-parallel_sjob` alone (which produces the
PDF on its own).

### Stale Sserver across sessions                           {#valid-run-gotcha-sserver}

`trustutils.run` auto-detects a running Sserver and uses it whether or not
you asked for one. A leftover Sserver from a previous run, possibly with
the wrong proc count, will silently affect subsequent reports. Always
stop the Sserver (`Sserver -1` or `Sserver -2`) at the end of a session.

### Form layout requirements                                 {#valid-run-gotcha-layout}

- `Run_fiche` must be invoked from the form root, never from `src/`.
- The form root must contain exactly one `.ipynb` file.

### `-post_run` is rejected for Jupyter forms                {#valid-run-gotcha-postrun}

`-post_run` is a PRM-only flag. For Jupyter forms, drive post-execution
work through Python `post_run` hooks on each case
(see @ref Valid_WriteReport).


## PRM conversion                                           {#valid-run-prm}

The historical PRM format is being discontinued. If you inherit a PRM
report and need to touch it, convert it to a Jupyter notebook first:

```
Run_fiche -convert_prm
```

This runs on a PRM report directory (the one containing `src/*.prm`) and
emits a Jupyter notebook alongside the original. Everything else in this
page then applies.


## Pointers                                                 {#valid-run-pointers}

- @ref Valid_WriteReport — authoring a new validation report.
- @ref Dev_ValidationFormAPI — doxygen API of `trustutils`.
- `$TRUST_ROOT/bin/Sjob/README` — Sserver reference.
- `Run_fiche -help` — full flag list (includes PRM-only flags not covered
  here).
