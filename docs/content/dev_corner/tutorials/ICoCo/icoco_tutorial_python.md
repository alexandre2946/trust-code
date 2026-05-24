@page DevTut_Icoco_python Using the ICoCo API from Python 

Since TRUST 1.8.3, the ICoCo API is accessible from Python through a SWIG
wrapper. This makes it possible to write supervision scripts entirely in Python
rather than in C++, with full access to the MEDCoupling Python API for field
manipulation and to `mpi4py` for multi-process coordination.

---

## Setting Up the Python Environment

A dedicated environment script provides the necessary Python paths and MPI
bindings. It must be sourced before launching any Python-based ICoCo supervisor:

```bash
source $TRUST_ROOT/env_for_python.sh
```

Verify the availability of the three core packages:

```python
import trusticoco as ti
import medcoupling as mc
from mpi4py import MPI
```

If all three imports succeed without error, the environment is correctly configured.

---

## Running an ICoCo Case from Python

The Python ICoCo API mirrors the C++ API: the same methods (`initialize()`,
`initializeTimeStep()`, `solveTimeStep()`, `validateTimeStep()`, `terminate()`,
etc.) are available on the Python problem object returned by `trusticoco`.

The data file preparation for a Python-supervised run is identical to a C++
supervised run: the TRUST data file must declare the ICoCo problem name and must
not contain a `Solve pb` instruction.

### Example Case

A fully worked example covering the conjugate heat transfer case from the previous
tutorial is available in the TRUST test database:

```bash
source $TRUST_ROOT/env_for_python.sh
trust -copy docond_VEF_3D_ICoCo_py
```

Examine the `prepare` script in the copied directory for the complete setup
procedure and then study the `main.py` file to see how the Python API is used to
drive two coupled problems.

@attention Even in sequential mode (single-process), a Python ICoCo supervisor
must always be launched through `mpirun`:

```bash
mpirun -np 1 python main.py
```

Launching `python main.py` directly will fail because the TRUST problem object
requires an MPI communicator to be initialized.

---

## Using ICoCo with Python in a Baltik Project

For users who develop their own TRUST-based application through the baltik build
system, a dedicated repository provides examples and build scripts for integrating
the Python ICoCo wrapper into a custom baltik:

```
https://github.com/cea-trust-platform/icoco-swig-baltik
```

This repository documents how to add the `trusticoco` Python module as a
dependency of an existing baltik project and how to structure Python supervisors
that interact with custom physics modules.

---

## Summary

The Python interface to ICoCo provides a convenient alternative to C++ supervision
scripts, particularly when the coupling logic involves complex data transformations
that benefit from the expressiveness of Python and the MEDCoupling Python API.
The main points to remember are:

The `trusticoco` Python module exposes the same interface as the C++ ICoCo API.
The MEDCoupling Python API (`medcoupling`) is used for field construction,
interpolation, and manipulation. Multi-process runs require `mpi4py` and must
always be launched with `mpirun`, regardless of the number of processes. Data file
preparation follows the same rules as for C++ supervisors: declare
`ICoCoProblemName` and remove `Solve pb`.