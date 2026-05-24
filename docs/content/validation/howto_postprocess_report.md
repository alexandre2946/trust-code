@page Valid_PostprocessReport Postprocessing and presenting results

This page covers the **postprocessing half** of a Jupyter validation form:
loading the artifacts produced by `run.runCases()`, plotting and tabulating
them, including VisIt renderings, and — most importantly — presenting the
results in a way that future maintainers can actually act on.

Test-case declaration (everything that runs *before* `run.runCases()`) is
covered in @ref Valid_WriteReport. Running a form is covered in
@ref Valid_RunReport.

[TOC]

## Where postprocessing fits in the form                    {#valid-post-where}

A validation form has a strict ordering:

```
markdown + run.introduction + run.TRUST_parameters
case declarations (run.addCase, addCaseFromTemplate, ...)
run.runCases()
─── everything below this line is postprocessing ───
plot/visit/table cells, narrative, displayed images
run.tablePerf()
```

The boundary is not just stylistic. The non-regression extraction process
concatenates the Python cells of the notebook **up to and including the
`run.runCases()` call**, replaces that call with a list-emitting helper, and
runs the result in a headless context. Anything that lives below that line
is invisible to NR; anything that lives above must be safe to execute in
that headless mode. Several of the gotchas below follow from that single
rule.


## Code guidelines                                          {#valid-post-code}

### All artifacts live under `run.BUILD_DIRECTORY`          {#valid-post-buildir}

Every file the postprocessing code reads or writes must resolve to a path
inside `run.BUILD_DIRECTORY`. A validation form can be redirected to a
different build directory (e.g. by passing `-dest` through
`JUPYTER_RUN_OPTIONS`), so hardcoded `./build` paths break in those modes.

Most of the `trustutils` API does the right thing automatically: relative
paths given to `plot.loadText`, `Graph.addPoint`, `Graph.addSegment`,
`Table.addLine`-with-a-file, `dumpText`, `visit.showField`, etc. are
resolved against `run.BUILD_DIRECTORY`. Use those entry points when
possible.

For raw Python — anything that opens or writes a file directly, the most
common offender being `plt.savefig` — join with `run.BUILD_DIRECTORY`
explicitly:

```python
import os, matplotlib.pyplot as plt
plt.savefig(os.path.join(run.BUILD_DIRECTORY, "my_fig.png"))
```

If a path needs to be validated explicitly, `run.sanitizePathToBUILD_DIRECTORY(path)`
returns it as an absolute path inside `build/` (and raises if the input
escapes it).


### Make artifacts archive-friendly                          {#valid-post-archive}

`Run_fiche -archive` packages `build/` plus a `used_files` manifest that
lists every file the report depends on (see @ref Valid_RunReport). Files
not in `used_files` are not in the archive, and the report cannot be
rebuilt from it.

The `trustutils` API marks the files it touches automatically: every
`plot.loadText`, `plot.read_csv`, `Graph.addPoint/addSegment`,
`plot.displayImageFile`, `dumpText`, ... calls `run.saveFileAccumulator`
under the hood. **Custom code** that writes or reads a file outside the
`trustutils` API must do that itself:

```python
out = os.path.join(run.BUILD_DIRECTORY, "summary.csv")
np.savetxt(out, data)
run.saveFileAccumulator("summary.csv")     # mark for archiving
```

A form is PR-ready only when `Run_fiche -archive` succeeds after a full
build (see @ref valid-run-archive).


### No `matplotlib` or `trustutils.plot` before `run.runCases()`  {#valid-post-mpl-order}

Never `import matplotlib.*` and never `from trustutils import plot` (which
transitively imports matplotlib) in a cell that runs **before** the
`run.runCases()` call. NR extraction executes those cells in a headless
context, and the matplotlib import there can break extraction.

Concretely: keep all `from trustutils import plot`, `import matplotlib.*`
and `from trustutils import visit` imports in **postprocessing cells only**,
strictly below `run.runCases()`. If a helper in `src/python_modules/` uses
matplotlib, it must also only be imported from postprocessing cells.


### Other gotchas                                            {#valid-post-misc}

- **Static images in markdown.** Use `![](src/fig.png)` (markdown) to make
  the image appear in the PDF. The HTML form `<img src="src/fig.png">`
  renders in JupyterLab but does **not** land in the PDF.
- **VisIt-disabled mode.** Guard VisIt-dependent code with
  `if not visit.VisItDisabled(): ...` so the report still builds under
  `Run_fiche -no_visit` and on machines where VisIt is not available.
- **`run.saveFormOutput()`.** Opt in to keeping cell outputs (plots,
  tables, printed text) in the saved `.ipynb` file. Not needed for
  `Run_fiche -export_pdf`, which re-executes the notebook anyway, but
  useful when sharing a pre-executed notebook.


## Presentation guidelines                                  {#valid-post-present}

A validation form is read more times than it is written. The first reader
is yourself. The second is your colleagues. The third — the one that
matters most for this section — is a maintainer some years from now who
may not share your numerical-methods or physics background.

The cost of getting this wrong differs by project. In TRUST itself, which
is essentially a numerical toolbox, a thin validation form is annoying but
survivable: another expert can usually recover the intent. In **baltiks
that implement specific physical models, correlations or numerical
schemes**, the validation form is often the **only living explanation** of
those choices, so a thin form is much more expensive — the knowledge
attached to it walks out the door with whoever wrote it. *TrioCFD is one
open-source example of such a project.* The guidelines below are written
with that harder case in mind; they apply to any validation form, but they
matter more the further you are from a pure numerical toolbox.

The rules below are deliberately opinionated and a bit abstract. Section
[A practical checklist](#valid-post-checklist) at the end of this section
translates them into something a reviewer can tick off.


### Tell the reader what is asserted                         {#valid-post-assert}

For every figure, table or printed result, the surrounding markdown must
state two things: **what claim is being made**, and **why the visual
supports it**. A figure with no commentary is decoration, not validation;
if you cannot articulate why a visual is in the form, remove it.

Add a one-sentence summary at the top of the postprocessing section that
names the claim being validated — e.g. *"this report validates the
second-order convergence of scheme X on a manufactured solution"* — so
that a reader who only skims the form still gets the headline. Do not
assume the reader will infer intent from a graph.


### Make pass/fail visible                                   {#valid-post-passfail}

The strongest form of validation is a quantitative check against an
**absolute reference** (analytical solution, manufactured solution,
conservation law, asymptotic convergence rate, ...). When you have one,
materialise the comparison as a number with a stated tolerance — in a
table or as an annotation on the graph — not just as overlapping curves.
"The curves look the same" is not a pass criterion: a maintainer
regenerating the PDF cannot tell whether a small deviation is noise or a
regression.

When no absolute reference is available — a common situation in CFD — the
section must instead **explicitly describe the expected behaviour and the
plausible failure modes**:

- the qualitative shape of the curve or field (monotonicity, sign,
  scale, asymptotic limits),
- the most likely failure modes you can think of (drift, oscillation,
  asymmetry, wrong slope, mass loss, ...),
- if known, the historical bugs that broke this case and how they
  manifested in the output.

The goal is that a maintainer regenerating the report after a change can
decide **on their own** whether the output still passes — without
contacting the original author.


### Write for someone who doesn't share your background      {#valid-post-audience}

A maintainer three years from now may have a different numerical methods
background, a different physics background, or both. Therefore:

- Describe the geometry and the physical setup briefly; do not assume
  "everyone knows this case".
- State units, sign conventions and non-dimensionalisations explicitly.
- Cite reference papers where relevant — and **summarise the relevant
  claim from each**. *"As discussed in [Smith 2018]"* without context is
  useless to someone who cannot quickly get the paper.
- Define acronyms and project-specific jargon the first time they appear,
  even the ones that feel obvious.


### Validation vs regression                                 {#valid-post-vs-regression}

A form that only shows numbers matching numbers is a regression test
wearing the clothes of a validation report. It will tell you when
something changes; it will not tell anyone whether the change is good or
bad.

A form that also explains what the numbers mean, why they should match,
and what a future failure would look like, is documentation that doubles
as a regression test. The latter is what carries a project across team
changes — especially in baltiks where the form is often the de facto
specification of a physical model or correlation.

This is not a checklist item; it is the question to ask yourself before
merging: *if I left the project tomorrow, would someone unfamiliar with
this test be able to maintain it from the report alone?* If not, the
report needs more text, not more plots.


### A practical checklist                                    {#valid-post-checklist}

Six items to run through before opening a PR with a new or modified
validation form:

1. **Headline.** A one-sentence summary at the top of the postprocessing
   section names what is being validated.
2. **Per-visual rationale.** Each figure, table or printed result is
   surrounded by text saying what it shows and what it asserts.
3. **Identifiable pass/fail.** Either a quantitative check against an
   absolute reference with a stated tolerance, or an explicit description
   of expected behaviour and plausible failure modes. Not just curve
   overlap.
4. **Self-contained context.** Geometry, physical setup, units, sign
   conventions and non-dimensionalisations are stated, not assumed.
5. **References summarised.** Cited papers are summarised in one or two
   sentences; acronyms are defined on first use.
6. **No decoration.** Anything that cannot justify its presence with one
   paragraph is removed.


## Basic API usage                                          {#valid-post-api}

A concise tour of the `trustutils` postprocessing entry points. All accept
paths relative to `run.BUILD_DIRECTORY` and call
`run.saveFileAccumulator` under the hood, so files they touch are
automatically archive-friendly.

### Plotting and tables — `trustutils.plot`                  {#valid-post-api-plot}

- `plot.Graph(title, ...)` — main plotting entry point; supports
  `addPoint(son_file)` (point probe), `addSegment(son_file)` (segment
  probe), `addResidu(dt_ev_file)` (residual curves), `add(x, y, ...)`
  (arbitrary curve), plus `label`, `legend`, `scale`, `visu`, `addPlot`
  for multi-panel layouts.
- `plot.Table([col1, col2, ...])` — tabular output; `addLine`, `sum`.
- `plot.loadText(file, ...)` — wraps `np.loadtxt`, relative to
  `BUILD_DIRECTORY`.
- `plot.read_csv(file, ...)` — wraps `pandas.read_csv`, same path
  convention.
- `plot.displayImageFile(path)` — display an image inline in the report
  **and** mark it as a build artifact. Use this for any image you produced
  yourself with `plt.savefig` (or similar) so that the image both appears
  in the PDF and survives `Run_fiche -archive`.

### VisIt renderings — `trustutils.visit`                    {#valid-post-api-visit}

- `visit.showMesh(lata_file, mesh)` — quick mesh rendering.
- `visit.showField(lata_file, plottype, field, ...)` — quick field
  rendering (`Pseudocolor`, `Contour`, `Vector`, ...).
- `visit.Show(...)` — composite renderings: multiple fields, multi-panel
  layouts, slices, custom VisIt commands.
- `visit.export_lata_base(...)` — extract numerical values (min, max,
  profiles) from a lata file for use in tables or further plots.
- `visit.VisItDisabled()` — guard for `-no_visit` and VisIt-less
  environments (see [Other gotchas](#valid-post-misc)).

### Interactive fine-tuning — `trustutils.widget`            {#valid-post-api-widget}

- `widget.interface(lata_file, plottype, field)` — interactive VisIt
  visualisation tweaker. Use **only** while authoring; replace with the
  equivalent static `visit.Show(...)` call once the visualisation is
  dialled in. Widgets in the final form are heavy and not portable.

### Helpers from `trustutils.run`                            {#valid-post-api-run}

- `run.BUILD_DIRECTORY` — absolute path of the build directory.
- `run.saveFileAccumulator(path)` — mark a file as a build artifact
  (called automatically by the entry points above; use it for custom
  files).
- `run.sanitizePathToBUILD_DIRECTORY(path)` — resolve and validate a path
  against `BUILD_DIRECTORY`.
- `run.extraitCoupe(data_file, probe_name, ...)` — extract probe values
  via the `extrait_coupe` tool.
- `run.dumpText(file, list_keywords=[...])` — print a text file inline in
  the report with optional keyword highlighting.
- `run.tablePerf()` — performance summary across all run cases; conventional
  last cell of a form.

For end-to-end examples of all of the above, see
`Validation/Rapports_automatiques/Verification/SampleFormJupyter`.


## Pointers                                                 {#valid-post-pointers}

- @ref Valid_WriteReport — authoring the test-case declaration half.
- @ref Valid_RunReport — running and exporting a report (`-archive` PR check).
- @ref Dev_ValidationFormAPI — doxygen API of `trustutils`.
- `Validation/Rapports_automatiques/Verification/SampleFormJupyter` —
  showcase form covering every API entry point above.
