@page DevTut_Icoco_basics Environment Setup and Basic Test Case 

This page covers the TRUST and ICoCo environment initialization, then walks through
a four-step exercise on a single-problem test case. The goal is to become familiar
with the ICoCo workflow before moving on to genuinely coupled problems.

---

## TRUST Environment Setup

Before starting, the TRUST environment must be sourced. On CEA Saclay workstations,
specific TRUST versions are made available via a dedicated script:

```bash
source /home/triou/env_TRUST_1.9.5.sh
```

On a personal workstation where TRUST has been installed locally, source the
environment from the installation root:

```bash
source $MyPathToTRUSTversion/env_TRUST.sh
```

To verify that the configuration is correct and to locate the TRUST root directory:

```bash
echo $TRUST_ROOT
```

## ICoCo Environment Setup

Source the ICoCo environment script bundled with TRUST:

```bash
source $TRUST_ROOT/Outils/ICoCo/ICoCo_src/full_env_MEDICoCo.sh
```

Check whether the ICoCo library has already been compiled:

```bash
ls $exec
```

If the command reports that the path does not exist (e.g., `ls: cannot access
.../ICoCo_opt: No such file or directory`), the library must be built first:

```bash
cd $TRUST_ROOT/Outils/ICoCo/ICoCo_src
baltik_build_configure -execute
make optim debug
source full_env_MEDICoCo.sh
```

@note If you do not have write access to the TRUST installation directory — typically
because TRUST was installed by a system administrator — copy the ICoCo source tree
to a local directory before building:

```bash
mkdir -p ICoCo
cp -r $TRUST_ROOT/Outils/ICoCo/ICoCo_src ICoCo/ICoCo_src
cd ICoCo/ICoCo_src
baltik_build_configure -execute
make optim debug
source full_env_MEDICoCo.sh
```

---

## Basic Test Case

This exercise uses the `Vahl_Davis_hexa` benchmark, a natural convection problem
discretized on a hexahedral mesh. The objectives are:

1. Run the case directly with the TRUST executable to establish a reference solution.
2. Adapt the case to run sequentially through ICoCo without any field exchange, and
   verify that results are bit-for-bit identical.
3. Introduce a first field exchange: impose a boundary condition value from the
   supervisor rather than from the data file.
4. Run the ICoCo-driven simulation in parallel on two MPI processes.

### Step 1 — Reference Simulation with TRUST

Create a working directory and copy the test case from the TRUST test database:

```bash
mkdir ICoCo_exercises
cd ICoCo_exercises
trust -copy Vahl_Davis_hexa
mv Vahl_Davis_hexa Vahl_Davis_hexa_trust
cd Vahl_Davis_hexa_trust
mv Vahl_Davis_hexa.data Vahl_Davis_hexa_trust.data
trust Vahl_Davis_hexa_trust
```

This produces a `Vahl_Davis_hexa_trust.lml` result file, which will serve as the
reference for all subsequent comparisons in this exercise.

### Step 2 — Sequential ICoCo Run Without Field Exchange

#### Adapting the Data File

Return to the exercises directory and create a separate copy of the case to be
driven through ICoCo:

```bash
cd ..
trust -copy Vahl_Davis_hexa
mv Vahl_Davis_hexa Vahl_Davis_hexa_ICoCo
cd Vahl_Davis_hexa_ICoCo
mv Vahl_Davis_hexa.data Vahl_Davis_hexa_ICoCo.data
```

Edit `Vahl_Davis_hexa_ICoCo.data` and apply two modifications. First, insert the
following line immediately after the `dimension` declaration to register the problem
with ICoCo:

```
Nom ICoCoProblemName Lire ICoCoProblemName pb
```

Second, comment out the `solve pb` instruction at the end of the file. The time
loop is now delegated to the supervisor.

#### Creating the Supervisor

The supervisor is a small C++ program that instantiates the problem object and
drives the time loop. Copy the provided template:

```bash
file="main_Vahl_Davis_hexa_ICoCo.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

Examine `main.cpp` to understand how the `Problem` object is constructed, how
`initializeTimeStep()` / `validateTimeStep()` are called, and how the time loop
is controlled. This sequential version uses a single MPI process.

#### Generating the Makefile and Compiling

```bash
cp $project_directory/share/bin/create_Makefile .
sh create_Makefile
make
```

The build step produces an executable named `couplage` together with a
`couplage.data` configuration file.

#### Running and Comparing Results

```bash
./couplage
compare_lata Vahl_Davis_hexa_ICoCo.lml \
             ../Vahl_Davis_hexa_trust/Vahl_Davis_hexa_trust.lml
```

The two result files must be identical: without any field exchange, the ICoCo
supervisor simply reproduces the TRUST standalone run.

### Step 3 — First Field Exchange: Supervisor-Imposed Boundary Condition

This step demonstrates how the supervisor can push values into the simulation at
runtime by writing to a named input field declared in the data file.

#### Setting Up the New Case

```bash
cd ..
cp -r Vahl_Davis_hexa_ICoCo Vahl_Davis_hexa_ICoCo_exchange
cd Vahl_Davis_hexa_ICoCo_exchange
trust -clean
mv Vahl_Davis_hexa_ICoCo.data Vahl_Davis_hexa_ICoCo_exchange.data
```

#### Adapting the Data File

In `Vahl_Davis_hexa_ICoCo_exchange.data`, locate the boundary condition on the
left wall and replace the static imposed temperature:

```
Gauche Paroi_temperature_imposee Champ_Front_Uniforme 1 10.
```

with an ICoCo input field declaration:

```
Gauche Paroi_temperature_imposee ch_front_input { nb_comp 1 nom TEMPERATURE_IN_DOM probleme pb }
```

The field name `TEMPERATURE_IN_DOM` must match the name used in the supervisor
when constructing and pushing the input field.

#### Updating the Supervisor

Copy the exchange-capable supervisor:

```bash
file="main_Vahl_Davis_hexa_ICoCo_exchange.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

#### Building, Running, and Comparing

```bash
make
./couplage
compare_lata Vahl_Davis_hexa_ICoCo_exchange.lml \
             ../Vahl_Davis_hexa_trust/Vahl_Davis_hexa_trust.lml
```

@note Results match the reference at all time steps except the first. This is
expected: the input field is written by the supervisor after the problem has been
initialized, so the very first time step applies the supervisor-provided value
rather than the default initial condition used in the reference run. This transient
discrepancy disappears from the second time step onward.

### Step 4 — Parallel Run on Two MPI Processes

#### Setting Up the Parallel Case

```bash
cd ..
cp -r Vahl_Davis_hexa_ICoCo_exchange Vahl_Davis_hexa_ICoCo_para
cd Vahl_Davis_hexa_ICoCo_para
make clean
rm main.cpp Vahl_Davis_hexa_ICoCo_exchange.lml
mv Vahl_Davis_hexa_ICoCo_exchange.data Vahl_Davis_hexa_ICoCo_para.data
```

Partition the domain and generate the parallel data files:

```bash
trust -partition Vahl_Davis_hexa_ICoCo_para
```

Then copy the parallel supervisor template:

```bash
file="main_Vahl_Davis_hexa_ICoCo_para.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

Open `main.cpp` and examine how processor subsets are assigned (search for
`dom_ids`) and how data file names are specified (search for `data_file`). In a
parallel ICoCo run, each participating problem is assigned a contiguous subset of
MPI ranks.

#### Building and Running

```bash
sh create_Makefile
make
mpirun -np 2 ./couplage
```

#### Comparing Results

```bash
compare_lata PAR_Vahl_Davis_hexa_ICoCo_para.lml \
             ../Vahl_Davis_hexa_ICoCo_exchange/Vahl_Davis_hexa_ICoCo_exchange.lml
```

The sequential and parallel runs must produce the same results, confirming that
the ICoCo supervisor correctly manages domain decomposition.