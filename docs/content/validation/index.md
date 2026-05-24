@page TRUST_Validation TRUST Validation

This section aims at explaining the validation process of TRUST.
It is directed at both users who are interested in how we validate the project,
and developers who need to write new validation tests.


## What "validation" means in TRUST                          {#valid-index-what}

TRUST's validation is **artifact-driven**. Each *validation form* covers a
specific physical or numerical case and demonstrates that the code
produces a correct result on it, compared against a documented reference
(analytical solution, established benchmark, experimental data, ...). A
form lives next to the source code, is written and maintained by
developers, and produces a PDF report that any user can read.

Validation is **not** a one-shot snapshot: the same forms are
re-executed automatically as non-regression tests, so a behaviour change
that breaks a previously-validated case is caught the moment it is
introduced. Reading the report of a form tells you what TRUST is known
to do well on a given case; the fact that the same form keeps running
tells you the behaviour is still current.


## One pipeline for TRUST and its baltiks                    {#valid-index-baltiks}

TRUST itself is a CFD kernel: schemes, discretizations, solvers,
parallelism. A **baltik** is an application code built on top of TRUST,
typically adding specific physical models, correlations or
application-level workflows. [TrioCFD](https://triocfd.cea.fr/) is one
open-source example of a baltik.

Both TRUST itself and every baltik ship their own validation forms,
using the same notebook-based pipeline, at different layers of
abstraction. For a user investigating a feature in a baltik, this means
the relevant validation evidence may be split across two layers: the
**TRUST-level forms** validate the numerical correctness of the kernel
(schemes, solvers, parallelism), and the **baltik-level forms** validate
the physical models or correlations layered on top. Both are usually
worth consulting.


## Where reports live                                        {#valid-index-where}

Validation forms live under `Validation/Rapports_automatiques/<Category>/<FormName>/`
at the project root — same convention in TRUST and in every baltik. Each
form ships its source (Jupyter notebook + dataset templates + meshes); no
pre-built PDFs are published. You generate the PDF report locally by
running `Run_fiche` in the form directory, which produces
`build/rapport.pdf`. The build command set and its options are documented
in @ref Valid_RunReport.

To find a form covering a topic, browse the category directory names
under `Rapports_automatiques/`; e.g. `Verification/` collects
correctness checks against analytical or manufactured references.


## What a report contains                                    {#valid-index-anatomy}

When you open a generated `rapport.pdf`, you will typically encounter,
in order:

- an **introduction** stating what is being validated and against which
  reference,
- a **test-case description**: geometry, physics, boundary and initial
  conditions, the numerical scheme(s) used,
- the **execution** of the test case(s),
- the **results**: curves, tables, field snapshots, along with the
  comparison against the reference,
- a short **conclusion** stating whether the agreement is acceptable,
- a final **performance summary** recording the cost of the runs.


## Detailed pages                                            {#valid-index-pages}

**For users and developers** — running an existing report from source:

- @subpage Valid_RunReport

**For developers** — authoring a new validation form:

- @subpage Valid_WriteReport — declaring the test cases.
- @subpage Valid_PostprocessReport — postprocessing and presenting results.
- @subpage Dev_ValidationFormAPI — full doxygen API of `trustutils`.
