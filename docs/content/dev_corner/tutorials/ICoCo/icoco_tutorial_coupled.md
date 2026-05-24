@page DevTut_Icoco_coupled Coupled Problem with ICoCo 

This exercise builds on the single-problem workflow established in the previous
tutorial and addresses the core use case of ICoCo: coupling two physically distinct
problems that exchange field data at runtime. The test case used is
`docond_VEF_3D`, a conjugate heat transfer problem in which a solid domain and a
fluid domain are thermally coupled through a shared interface.

The exercise is structured in four steps:

1. Run the fully coupled reference case directly with TRUST to establish a reference
   solution.
2. Split the monolithic data file into two separate problem files and drive them
   together through ICoCo without any field exchange.
3. Introduce a **one-way coupling**: the fluid problem computes its temperature field
   independently, and the resulting boundary temperature is imposed on the solid
   problem at each time step.
4. Introduce a **two-way coupling**: the solid problem additionally sends its surface
   heat flux back to the fluid problem, implementing the physically consistent
   conjugate heat transfer algorithm.

---

## Step 1 — Reference Coupled Simulation with TRUST

### Sequential and Parallel Reference Runs

Copy the coupled reference case from the TRUST test database and run it both
sequentially and in parallel:

```bash
cd ICoCo_exercises
trust -copy docond_VEF_3D
mv docond_VEF_3D docond_VEF_3D_trust
cd docond_VEF_3D_trust
trust docond_VEF_3D
trust -partition docond_VEF_3D
trust PAR_docond_VEF_3D 2
```

Verify that the sequential and parallel results are consistent:

```bash
compare_lata docond_VEF_3D.lml PAR_docond_VEF_3D.lml
```

Differences between sequential and parallel results should be well below the
tolerance threshold of 10⁻⁵, as expected for a correctly partitioned computation.

---

## Step 2 — Splitting the Problem for ICoCo

The monolithic `docond_VEF_3D` data file must be decomposed into separate files:
one defining the mesh geometry and partitioning, and two independent problem files
(one per domain) that will be driven by the ICoCo supervisor.

### Creating the ICoCo Case Directory

```bash
cd ..
trust -copy docond_VEF_3D
mv docond_VEF_3D docond_VEF_3D_ICoCo
cd docond_VEF_3D_ICoCo
```

### Generating Partitioned Meshes

Create a mesh data file for the solid domain by copying the full data file and
retaining only the geometry and partitioning blocks up to (and including) the
`End` instruction. Remove the time scheme, the problem definitions, and any
post-processing blocks. Uncomment the partition block:

```bash
cp docond_VEF_3D.data docond_VEF_3D_mesh1.data
```

Edit `docond_VEF_3D_mesh1.data` so that it contains only the geometry, mesh, and
partitioning information for the solid domain (`DOM1`). Similarly, create a mesh
file for the fluid domain:

```bash
cp docond_VEF_3D_mesh1.data docond_VEF_3D_mesh2.data
```

Edit `docond_VEF_3D_mesh2.data` so that it retains only the fluid domain (`DOM2`)
information. Run both mesh files to generate the partitioned zone files:

```bash
trust docond_VEF_3D_mesh1
trust docond_VEF_3D_mesh2
```

This produces four `.Zones` files that describe the domain decomposition:

```
DOM1_0000.Zones   DOM1_0001.Zones
DOM2_0000.Zones   DOM2_0001.Zones
```

### Creating Individual Problem Data Files

Create a data file for the solid domain's problem:

```bash
cp docond_VEF_3D.data docond_VEF_3D_dom1.data
```

In `docond_VEF_3D_dom1.data`:
- Remove the `Probleme_Couple pbc`, `Associate pbc pb1`, `Associate pbc pb2`
  lines and any per-problem file output specifiers.
- Replace `Associate pbc sch` with `Associate pb sch`.
- Replace `Discretize pbc dis` with `Discretize pb dis`.
- Replace `Solve pbc` with `Solve pb`.
- Remove the mesh and partitioning blocks, which are now handled by the mesh
  files, and uncomment the `scatter` block.
- In the data file, retain only the information related to the solid domain
  (`pb1`, `dom_solide`) and perform a global substitution `pb1` → `pb`.

Create the fluid domain's problem file analogously:

```bash
cp docond_VEF_3D_dom1.data docond_VEF_3D_dom2.data
```

In `docond_VEF_3D_dom2.data`, retain only the fluid domain information
(`pb2`, `dom_fluide`) and substitute `pb2` → `pb` throughout.

The two domains exchange heat through their shared interface: the `Paroi_echange1`
boundary in the solid domain and `Paroi_echange2` in the fluid domain. For this
uncoupled baseline step, replace the contact boundary conditions with a simple
fixed-temperature condition on both sides:

In `docond_VEF_3D_dom1.data`:
```
Paroi_echange1 paroi_temperature_imposee Champ_Front_Uniforme 1 50.
```

In `docond_VEF_3D_dom2.data`:
```
Paroi_echange2 paroi_temperature_imposee Champ_Front_Uniforme 1 50.
```

### Verifying Independent Problem Runs

Before activating ICoCo, verify that each problem file runs correctly in isolation:

```bash
trust docond_VEF_3D_dom1 2
trust docond_VEF_3D_dom2 2
```

Both problems should complete successfully. At this stage there is no coupling:
each problem applies its own fixed boundary temperature independently.

Produce per-problem result files that will be used later for comparison. In the
`Post_processing` block of each problem, specify distinct output file names:
- add `fichier pb1` in the solid problem's post-processing block,
- add `fichier pb2` in the fluid problem's post-processing block.

Rerun the coupled TRUST case to regenerate these labelled output files:

```bash
trust docond_VEF_3D 2
```

### Adapting the Data Files for ICoCo

Add the ICoCo problem registration line immediately after `dimension 3` in both
`docond_VEF_3D_dom1.data` and `docond_VEF_3D_dom2.data`:

```
Nom ICoCoProblemName Lire ICoCoProblemName pb
```

Remove the `Solve pb` instruction from both files; the time advancement is now
delegated to the ICoCo supervisor.

### Running with ICoCo (No Exchange)

Copy the provided supervisor for this step:

```bash
file="main_docond_VEF_3D_ICoCo.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

Open `main.cpp` and identify:
- where MPI processor subsets are declared (`dom_ids`),
- where the data file names are specified (`data_file`),
- the time-stepping loop (`while`).

Build and run on four MPI processes:

```bash
sh $project_directory/share/bin/create_Makefile 4
make
mpirun -np 4 ./couplage
```

Compare the ICoCo results with the independent TRUST runs:

```bash
compare_lata pb1.lml docond_VEF_3D_dom1.lml
compare_lata pb2.lml docond_VEF_3D_dom2.lml
```

@note Differences are expected at this stage, because the contact boundary
condition from the reference TRUST run has been replaced by a fixed temperature of
50 °C on both sides. The purpose of this step is only to confirm that the ICoCo
supervisor correctly initializes and advances both problems simultaneously.

---

## Step 3 — One-Way Coupling

One-way coupling means that one problem (here the fluid) runs independently, and
the field it produces at the coupling interface is imposed as a boundary condition
on the other problem (here the solid) at each time step. The fluid field drives the
solid, but not vice versa.

### Setting Up the Case

```bash
cd ..
cp -r docond_VEF_3D_ICoCo docond_VEF_3D_ICoCo_coupling1
cd docond_VEF_3D_ICoCo_coupling1
```

### Defining the Output Field in the Fluid Problem

In `docond_VEF_3D_dom2.data`, add the following field definition inside the
`Definition_champs` block of the post-processing section. This field extracts the
temperature at the fluid side of the coupling interface and makes it available for
export via ICoCo:

```
TEMPERATURE_OUT_DOM2 Interpolation {
    localisation elem
    domaine dom_fluide_boundaries_Paroi_echange2
    source refChamp { Pb_Champ pb temperature }
}
```

The domain name `dom_fluide_boundaries_Paroi_echange2` is an automatically
constructed name identifying the boundary `Paroi_echange2` of the fluid domain.

### Defining the Input Field in the Solid Problem

In `docond_VEF_3D_dom1.data`, replace the fixed boundary condition on
`Paroi_echange1` with an ICoCo input field:

```
Paroi_echange1 paroi_temperature_imposee ch_front_input {
    nb_comp 1 nom TEMPERATURE_IN_DOM1 probleme pb }
```

The field name `TEMPERATURE_IN_DOM1` must match the name used in the supervisor
when writing the interpolated fluid temperature into the solid boundary.

### Updating the Supervisor

Copy the one-way coupling supervisor:

```bash
file="main_docond_VEF_3D_ICoCo_coupling1.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

Compare this supervisor against the no-exchange version to understand the
additions:

```bash
diff main.cpp ../docond_VEF_3D_ICoCo/main.cpp
```

The key additions are a `TrioDEC` object that handles MPI data redistribution
between the processor subsets of the two problems, and a `TrioField` object that
carries the field values. The coupling sequence within the time loop is:

1. Solve the fluid problem for the current time step.
2. Retrieve `TEMPERATURE_OUT_DOM2` from the fluid problem via ICoCo.
3. Interpolate and redistribute the field data using the `TrioDEC`.
4. Push the result into the solid problem as `TEMPERATURE_IN_DOM1` via ICoCo.
5. Solve the solid problem for the current time step.

### Building, Running, and Comparing

```bash
make
mpirun -np 4 ./couplage
```

Compare results with the uncoupled ICoCo run:

```bash
compare_lata docond_VEF_3D_dom1.lml ../docond_VEF_3D_ICoCo/docond_VEF_3D_dom1.lml
compare_lata docond_VEF_3D_dom2.lml ../docond_VEF_3D_ICoCo/docond_VEF_3D_dom2.lml
```

The fluid domain (`dom2`) results are unchanged, since the fluid problem runs
independently. The solid domain (`dom1`) results differ, confirming that the
coupling is active and that the solid boundary condition is now driven by the fluid
temperature.

---

## Step 4 — Two-Way Coupling

Two-way coupling for conjugate heat transfer requires using complementary boundary
condition types on each side of the interface: a Dirichlet condition (imposed
temperature) on one side and a Neumann condition (imposed heat flux) on the other.
Imposing Dirichlet conditions on both sides simultaneously is mathematically
ill-posed and must be avoided.

In this step, the solid problem continues to receive the fluid temperature as a
Dirichlet condition (as in Step 3), while the fluid problem now receives the
surface heat flux from the solid as a Neumann condition.

@warning Using Dirichlet boundary conditions on both sides of a thermal coupling
interface leads to an ill-posed problem and will generally produce non-physical
results or convergence failure. Always combine Dirichlet on one side with Neumann
on the other.

### Setting Up the Case

```bash
cd ..
cp -r docond_VEF_3D_ICoCo_coupling1 docond_VEF_3D_ICoCo_coupling2
cd docond_VEF_3D_ICoCo_coupling2
```

### Defining the Output Heat Flux Field in the Solid Problem

In `docond_VEF_3D_dom1.data`, add the following field definition in the
`Definition_champs` block of the post-processing section. It extracts the
surface heat flux at `Paroi_echange1` and exposes it via ICoCo:

```
FLUX_SURFACIQUE_OUT_DOM1 Interpolation {
    localisation elem
    domaine dom_solide_boundaries_Paroi_echange1
    source Morceau_equation {
        type operateur numero 0 option flux_surfacique_bords
        source refChamp { Pb_Champ pb temperature }
    }
}
```

The domain name `dom_solide_boundaries_Paroi_echange1` is the automatically
constructed name for the boundary `Paroi_echange1` of the solid domain.

### Defining the Input Heat Flux Field in the Fluid Problem

In `docond_VEF_3D_dom2.data`, replace the fixed boundary condition on
`Paroi_echange2` with an ICoCo Neumann input field:

```
Paroi_echange2 paroi_flux_impose ch_front_input {
    nb_comp 1 nom FLUX_SURFACIQUE_IN_DOM2 probleme pb }
```

### Updating the Supervisor

Copy the two-way coupling supervisor or start from the one-way version and add
the reverse exchange:

```bash
file="main_docond_VEF_3D_ICoCo_coupling2.cpp"
cp $TRUST_ROOT/doc/TRUST/exercices/ICoCo/$file main.cpp
```

Compare against the one-way supervisor to identify the additions:

```bash
diff main.cpp ../docond_VEF_3D_ICoCo_coupling1/main.cpp
```

The additions are a second `TrioDEC` object (`dec_flux`) handling data transfer
from the solid processor subset to the fluid subset, and a corresponding
`TrioField` (`field_flux`). The time loop now implements the following exchange
sequence:

1. Solve the fluid problem.
2. Export `TEMPERATURE_OUT_DOM2` from the fluid and import it as
   `TEMPERATURE_IN_DOM1` into the solid (as in Step 3).
3. Solve the solid problem.
4. Export `FLUX_SURFACIQUE_OUT_DOM1` from the solid and import it as
   `FLUX_SURFACIQUE_IN_DOM2` into the fluid.

### Building and Running

```bash
make
mpirun -np 4 ./couplage
```

### Discussion

Compare the two-way coupling results against the reference TRUST coupled run:

```bash
compare_lata docond_VEF_3D_dom1.lml ../docond_VEF_3D_ICoCo/pb1.lml
compare_lata docond_VEF_3D_dom2.lml ../docond_VEF_3D_ICoCo/docond_VEF_3D_dom2.lml
```

Transient differences are expected. The TRUST internal coupling uses an implicit
treatment of the interface condition, whereas the ICoCo-based coupling applies an
explicit staggered scheme: values exchanged at time step *n* are used to compute
the solution at time step *n+1*. As the simulation progresses toward a steady
state, the two approaches converge toward the same solution.

@note For a detailed analysis of the convergence behaviour and a comparison between
implicit and explicit conjugate heat transfer coupling strategies, refer to the
TRUST validation report `CouplageFluideSolide`, available within the installation
at `$project_directory/share/Validation/Rapports_automatiques/CouplageFluideSolide`.