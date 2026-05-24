@page GG_Intro Introduction

**TRUST** is a High Performance Computing (HPC) thermal-hydraulic engine for Computational Fluid Dynamics (CFD) developed at the Department of %System and Structure Modelling (DM2S) of the French Atomic Energy Commission (CEA).

The acronym **TRUST** stands for <b>TR</b>io_<b>U</b> <b>S</b>oftware for <b>T</b>hermalhydraulics. This software was originally designed for conduction, incompressible single-phase, and Low Mach Number (LMN) flows with a robust Weakly-Compressible (WC) multi-species solver. However, a huge effort has been conducted recently, and now TRUST is able to simulate real compressible multi-phase flows.

TRUST is also being progressively ported to support GPU acceleration (NVidia/AMD), using the [Kokkos](https://kokkos.org/kokkos-core-wiki/) library.

The software is open source with a [BSD license](https://github.com/cea-trust-platform/trust-code/blob/master/License.txt), available on GitHub via [this link](https://github.com/cea-trust-platform/trust-code).

You can easily create new projects based on the **TRUST** platform. These projects are named **BALTIK** projects (<b>B</b>uild an <b>A</b>pplication <b>L</b>inked to <b>T</b>RUST <b>K</b>ernel).

## A bit of history: the Modular Software Named Trio_U

**TRUST** was born from the splitting of the **Trio_U** software into two pieces. **Trio_U** was a software brick based on the **Kernel** brick (which contains the equations, space discretizations, numerical schemes, parallelism, ...) and was used by other CEA applications (see Figure 1).

![Figure 1: Trio_U brick software](figures/tikz1.png)

In 2015, **Trio_U** was divided into two parts: **TRUST** and **TrioCFD**.

- **TRUST** is a new platform; its name means: <b>TR</b>io_<b>U</b> <b>S</b>oftware for <b>T</b>hermohydraulics.
- **TrioCFD** is an open source BALTIK project based on **TRUST**.

Here are some other selected BALTIK projects based on the TRUST platform (see Figure 2).

![Figure 2: Selected BALTIK projects based on the TRUST platform.](figures/tikz2.png)

## Short History

**TRUST** is developed at the Laboratory of High Performance Computing and Numerical Analysis (LCAN) of the Software Engineering and Simulation Service (SGLS) in the Department of System and Structure Modeling (DM2S).

The project started in 1994 and improved versions have been built ever since:

- **1994:** Start of the project Trio_U
- **1997:** v1.0 - Finite Difference Volume (VDF) method only
- **1998:** v1.1 - Finite Element Volume (VEF) method only
- **2000:** v1.2 - Parallel MPI version
- **2001:** v1.3 - Radiation model (TrioCFD now)
- **2002:** v1.4 - LES turbulence models (TrioCFD now)
- **2006:** v1.5 - VDF/VEF Front Tracking method (TrioCFD now)
- **2009:** v1.6 - Data structure revamped
- **2015:** v1.7 - Separation TRUST & TrioCFD + switch to open source
- **2019:** v1.8 - New polyhedral discretization (PolyMAC) + Multiphase problem + Weakly Compressible model
- **2022:** v1.9.0 - Modern C++ code (templates, CRTP, ...) + remove MACROS + support GPU (NVidia/AMD)
- **2025:** v1.9.6 - Unified version to handle 32-64b integers + VEF discretization supported on GPU

## Where to go next

The rest of the General Guide walks through the practical steps of
writing and running a simulation:

- @ref GG_Data — anatomy of a `.data` file: blocks, objects, interpretors, and how to launch it
- @ref GG_Problems — choosing a problem type (hydraulic, thermohydraulic, multiphase, ...)
- @ref GG_BoundaryConditions — boundary conditions
- @ref GG_Solvers — pressure solvers
- @ref GG_TemporalSchemes — time integration schemes
- @ref GG_SpatialSchemes — spatial discretization schemes
- @ref GG_PostProcessing — output files and post-processing
- @ref GG_EndData — `Solve`, stop file, save and restart
- @ref GG_Parallel — partitioning and parallel execution

For step-by-step worked examples, see the @ref UserTutorials.

## TRUST command-line reference {#trust_cli_reference}

The `trust` script wraps the actual TRUST executable with a collection of
convenience options.  Run `trust -h` at any time for the live list — the
groups below are a quick overview.

<h3>Information</h3>

| Flag | What it does |
|---|---|
| `-help`, `-h` | List every option. |
| `-index` | Open the TRUST resource index. |

<h3>Project scaffolding</h3>

| Flag | What it does |
|---|---|
| `-baltik [name]` | Instantiate an empty baltik project. |
| `-eclipse-trust` | Generate Eclipse config to import TRUST sources. |
| `-eclipse-baltik` | Generate Eclipse config to import baltik sources (requires TRUST to have been configured under Eclipse first). |

<h3>Data-file authoring</h3>

| Flag | What it does |
|---|---|
| `-config nedit\|vim\|emacs\|gedit` | Set up your editor to syntax-highlight TRUST keywords. |
| `-edit` | Edit the data file. |
| `-jupyter` | Create a basic Jupyter notebook file. |
| `-copy` | Copy a test-case data file from the TRUST database into the current directory. |
| `-search KEYWORDS` | List test cases from the database that contain the given keywords. |

<h3>Mesh & parallel</h3>

| Flag | What it does |
|---|---|
| `-partition` | Partition the mesh (creates the `.Zones` files). |
| `-mesh` | Visualise the mesh(es) referenced by the data file. |
| `-c NCPUS` | Use NCPUS CPUs per task for a parallel calculation. |

<h3>Run monitoring & cleanup</h3>

| Flag | What it does |
|---|---|
| `-evol` | Monitor the calculation in a GUI. |
| `-probes` | Monitor probes only. |
| `-clean` | Remove every file TRUST generated in the current directory. |

<h3>Non-regression testing</h3>

| Flag | What it does |
|---|---|
| `-check all\|TESTCASE\|LISTFILE` | Run non-regression checks. |
| `-ctest all\|TESTCASE\|LISTFILE` | Same, via the CTest harness. |
| `-check function\|class\|class::method` | Check the tests covering a given function, class, or method. |
| `-ctest function\|class\|class::method` | Same, via CTest. |

<h3>Debugging & profiling</h3>

| Flag | What it does |
|---|---|
| `-gdb` | Run under gdb. |
| `-valgrind` | Run under valgrind. |
| `-valgrind_strict` | Run under valgrind with no suppressions. |
| `-callgrind` | Run the callgrind profiler (via valgrind). |
| `-massif` | Run the massif heap profiler (via valgrind). |
| `-heaptrack` | Run heaptrack (heap profile; faster than massif). |
| `-advisor` | Run Intel Advisor (vectorisation analysis). |
| `-vtune` | Run Intel VTune (profiling). |
| `-perf` | Run Linux `perf` (profiling). |
| `-trace` | Run the Intel Trace Analyzer (MPI profiling). |

<h3>Cluster submission</h3>

| Flag | What it does |
|---|---|
| `-create_sub_file` | Just create a submission file, don't submit. |
| `-prod` | Create and submit on the main production queue with exclusive resource. |
| `-bigmem` | Create and submit on the big-memory production queue. |
| `-queue QUEUE` | Create and submit to a specific queue. |

<h3>Passing options through to the underlying executable</h3>

| Flag | What it does |
|---|---|
| `datafile -help_trust` | Print the options of the underlying `TRUST_EXECUTABLE [CASE[.data]] [options]`. |
