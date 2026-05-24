@page DevTut_ValidationForm Validation form and non-regression test cases

A *validation form* is a Jupyter notebook that documents a test case —
its data file, mesh, expected results, plots and reference values.
TRUST baltiks ship them under
`share/Validation/Rapports_automatiques/`.  This tutorial walks through
two things you will do regularly while developing a baltik:

1. create a new validation form from an existing `.data` file;
2. promote that form into a non-regression test case that runs in your
   nightly checks.

The example uses a fictitious `upwind.data` case but the workflow is
identical for any test case you want to document.

## Create a new validation form

Each form lives in its own directory under
`share/Validation/Rapports_automatiques/`.  Make a fresh one and drop
the data file (and any companion files, such as the mesh) into a `src/`
subdirectory:

```bash
cd $project_directory
cd share/Validation/Rapports_automatiques
mkdir -p upwind/src

cp /path/to/upwind.data upwind/src/
cp /path/to/upwind.geo  upwind/src/
```

Generate the notebook skeleton with the `trust` script:

```bash
cd upwind
trust -jupyter
```

You now have an `upwind.ipynb` file in the directory — a stub Jupyter
notebook with the cells TRUST expects (data-file extraction, run,
post-processing).  Edit it to add the prose, plots and validation
checks that describe what this case is meant to demonstrate.

Add and commit:

```bash
git add upwind
git commit -m "New validation notebook"
```

## Run the form and build the report

Two helpers ship with the baltik environment:

```bash
cd $project_directory/share/Validation/Rapports_automatiques/upwind

Run_fiche                     # execute the notebook end-to-end
Run_fiche -export_pdf         # same, plus produce a PDF report
evince build/rapport.pdf &    # view it
```

The PDF will contain every field and probe post-processed by the data
file, plus the prose, plots and tables you wrote into the notebook.

## Promote the form to a non-regression test case

A validation form is a documentation artefact.  To make it run as part
of the baltik's nightly non-regression checks, you must register it in
`tests/Reference/Validation/`.  The `check_optim` target does this
automatically:

```bash
cd $project_directory
make check_optim
```

You should see lines like:

```
Creation of upwind_jdd1
Creation of upwind_jdd1/lien_fiche_validation
Extracting test case (upwind.data) ...End.
Creation of the file upwind_jdd1.lml.gz
```

A new test case `upwind_jdd1` is created under
`tests/Reference/Validation/`, pointing back at the validation form.

@note It is common for the **parallel** variant of a freshly-extracted
case (`PAR_<name>_jdd1`) to crash on first run with a `CORE` message.
The cause is almost always a `Scatter` path that points outside the
test-case directory (e.g. `Scatter ../upwind/DOM.Zones dom`).  See the
next section for the fix.

## Fixing a failing parallel test

When `make check_optim` reports a crash on `PAR_upwind_jdd1`, copy the
extracted case into the build dir, partition it manually, and
re-run:

```bash
cd $project_directory/build
trust -copy upwind_jdd1
```

If `trust -copy` complains that it doesn't know the case yet, re-run
`./configure` from the project root so the new test case is registered,
then retry the copy:

```bash
cd $project_directory
./configure
cd build
trust -copy upwind_jdd1
```

Now partition the case and run it on two processes to reproduce the
crash:

```bash
cd upwind_jdd1
trust -partition upwind_jdd1
trust PAR_upwind_jdd1 2
```

Edit `PAR_upwind_jdd1.data` to fix the `Scatter` line.  The common
correction is:

```text
- Scatter ../upwind/DOM.Zones dom
+ Scatter DOM.Zones dom
```

When the parallel run completes, propagate the fix back into the source
notebook directory so future extractions stay correct:

```text
$project_directory/share/Validation/Rapports_automatiques/upwind/src/upwind.data
```

Re-run only the previously-failing tests:

```bash
cd $project_directory
make check_last_pb_optim
```

You should see:

```
Successful tests cases: 1/1
```

## Commit the new test case

```bash
git status -uno
git add tests/Reference/Validation/upwind_jdd1/upwind_jdd1.data
git commit -m "New reference test"
```

## Run the full non-regression suite

Two convenience targets run every registered test:

```bash
make check_all_optim          # optimised binary
make check_all_debug          # debug binary
```

## Package the baltik for sharing

To produce a tarball that bundles the baltik (sources, validation
forms, registered test cases, build metadata):

```bash
make distrib
ls
```

The resulting `<baltik>.tar.gz` is what you hand to someone else who
wants to rebuild your baltik from scratch.
