@page Disc_PolyMAC_MPFA PolyMAC_MPFA

Unlike PolyMAC_CDO, PolyMAC_MPFA does not introduce the vorticity and does not require an explicit complex dual mesh. It is built on the Multi-Point Flux Approximation (MPFA) family of methods. The location of the unknowns is shown below.

![Location of the unknowns when using PolyMAC_MPFA](figures/PolyMAC_unknowns.png)

In the TRUST class hierarchy, `Domaine_PolyMAC_MPFA` and its operators inherit from their `PolyMAC_HFV` counterparts: the **only** ingredient that genuinely uses an MPFA reconstruction is the face-flux operator built by `Domaine_PolyMAC_MPFA::fgrad`. The convection operator, the mass solver and most of the Navier-Stokes machinery are inherited unchanged from PolyMAC_HFV.

## Degrees of freedom

| Field                              | Location                                                       |
|------------------------------------|----------------------------------------------------------------|
| Scalar (T, p, α, ...) — `Champ_Elem` | One value per cell                                             |
| Velocity — `Champ_Face`            | Normal component at each face + auxiliary tangential vector at each cell |

The cell-tangential auxiliary DOFs are *reconstructed* from the face-normal values after every linear solve (see [Cell-tangential velocity reconstruction](#cell-tangential-velocity-reconstruction) below). They are not independent unknowns.

## MPFA flux reconstruction

Three MPFA schemes are implemented:

- **MPFA-O** @cite AM08 @cite D14 — strongly consistent (exact on affine fields).
- **MPFA-O(η)** @cite ER98 — strongly consistent. The trace points are not at the face barycenter but are placed so as to make the local bilinear form symmetric, by solving a least-squares system with **LAPACK** `DGELSY`.
- **MPFA-symm** @cite LP05a @cite LP05b @cite LePotier2017 — symmetric and unconditionally coercive, but only **weakly** consistent: strongly consistent on simplices and Cartesian cells, and recovers consistency on general meshes by refinement arguments.

### Construction of the dual mesh

To build the sub-cell geometry on which the gradient is reconstructed, the procedure is:

1. Link each cell's barycenter (purple) to the barycenter of each cell's face \f$f\f$ (blue). This cuts each face into two sub-faces \f$\hat{f}_1\f$ and \f$\hat{f}_2\f$, subdividing each cell into \f$N_i\f$ sub-volumes \f$(S_{e,i})_{i\in\{1,\dots,N_i\}}\f$.
2. Introduce auxiliary trace unknowns \f$u_{f_v}\f$ at each sub-face (green). In MPFA-O they sit at the face barycenter; in MPFA-O(η) their offsets are chosen to symmetrise the local bilinear form (clamped to one-third / two-thirds of each half-edge).

![Construction of a gradient using the MPFA method](figures/MPFA.png)

On sub-volume \f$S_{e,i}\f$ the gradient of a potential \f$p\f$ is

\f[
G_{S_{e,i}}([p]_e) = \frac{1}{|S_{e,i}|}\sum_{f' \subset \partial S_{e,i}} (u_{f'_s} - p_e)\,\vec{n}_{f'}
\f]

with \f$\vec{n}_{f'}\f$ outward unit normals, assembled cell-by-cell as

\f[
G^{\text{MPFA}}_{|S_{e,i}} = G_{S_{e,i}}([p]_e),\qquad \forall e \in E,\; i \in S_e.
\f]

Continuity of the normal flux across each interior sub-face yields a small dense linear system relating the trace unknowns \f$u_{f_v}\f$ to the cell unknowns \f$u_e\f$. This system is solved (LAPACK `DGELSY`) and the trace unknowns are eliminated, leaving for each face \f$f\f$ a multi-point coefficient array `phif_c(f, e_s)` such that

\f[
[\vec{n}_f \cdot \nu \nabla u]_f \;\approx\; \sum_{e_s \in \mathrm{stencil}(f)} \texttt{phif\_c}(f, e_s)\,u_{e_s}.
\f]

### Per-vertex scheme selection

In TRUST the loop in `fgrad` is **over mesh vertices**: at each vertex, the small dense system above is built using the fan of incident cells and faces, and the result is accumulated into the global face stencil. For each vertex the three schemes are tried in order (MPFA-O → MPFA-O(η) → MPFA-symm); a candidate is **accepted** as soon as the smallest eigenvalue of the local symmetric bilinear form (`DSYEV`) satisfies

\f[
\lambda_{\min} > -10^{-8}\,|v|
\f]

where \f$|v|\f$ is the volume of the vertex sub-region. MPFA-symm is therefore invoked only as a **fallback** at vertices where both MPFA-O and MPFA-O(η) fail the coercivity test — typically in locally distorted or anisotropic regions of the mesh. At the first call of `fgrad`, the percentages of vertices using each scheme are printed, e.g.:

```
::fgrad(): 92.4% MPFA-O  6.3% MPFA-O(h)  1.3% MPFA-SYM
```

A high MPFA-SYM share does not necessarily indicate a numerical problem, but it tells the user where to look if the local convergence rate degrades: those are the vertices at which strong consistency had to be abandoned for guaranteed coercivity.

### Pressure / velocity reuse

`fgrad` carries an `is_p` flag that swaps the role of Dirichlet and Neumann boundary conditions inside the local solve, so the same routine reconstructs \f$[\vec{n}\cdot\nabla p]_f\f$ for the pressure-velocity coupling and \f$[\vec{n}\cdot\nu\nabla u]_f\f$ for diffusion without duplicate code.

## Diffusion operator

`Op_Diff_PolyMAC_MPFA_Elem` discretises

\f[
|e|\,\varphi_e\,\partial_t u_e \;-\; \sum_{f \in F_e} \sigma_{e,f}\,|f|\,[\vec{n}_f \cdot \nu \nabla u]_f \;=\; |e|\,s_e,
\f]

where \f$\sigma_{e,f} = \pm 1\f$ depending on face orientation. The face flux is the multi-point combination delivered by `fgrad`. On a K-orthogonal mesh the stencil collapses to the two-point TPFA pattern and recovers the classical orthogonal-mesh FV scheme; the percentage of faces where this collapse happens is printed by `Op_Diff_PolyMAC_MPFA_Elem::dimensionner_blocs` as a useful diagnostic of effective mesh orthogonality.

The same machinery is used by `Op_Diff_PolyMAC_MPFA_Face` for the viscous block of the momentum equation, applied componentwise to each scalar component of velocity.

## Incompressible Navier-Stokes

\f[
\begin{aligned}
\rho\varphi\,\partial_t u + \rho\,\nabla\cdot(u\otimes u) + \nabla p - \nabla\cdot\bigl(\nu(\nabla u + (\nabla u)^\top)\bigr) &= \rho f, \\
\nabla\cdot u &= 0.
\end{aligned}
\f]

**Mass equation** — discretised at the cell using Green-Ostrogradski:

\f[
|e|\,[\nabla\cdot u]_e = \sum_{f \in F_e} \sigma_{e,f}\,|f|\,[u]_f.
\f]

The face-normal flux depends only on face-normal velocities, which is the structural reason PolyMAC inherits the MAC-scheme div-free property on Cartesian meshes.

**Convective term** — `Op_Conv_EF_Stab_PolyMAC_MPFA_Face` is **face-based**: at each face \f$f\f$, the contribution of cell \f$e\f$ is

\f[
\frac{|f|\,u_f\,\varphi_e}{2}\,\bigl(1 + \alpha\,\mathrm{sgn}(u_f)\bigr)\,u_e,\qquad \alpha \in [0, 1],
\f]

where \f$\alpha\f$ blends pure upwind (\f$\alpha = 1\f$, exposed as `Op_Conv_Amont_PolyMAC_MPFA_Face`) and pure centred (\f$\alpha = 0\f$, exposed as `Op_Conv_Centre_PolyMAC_MPFA_Face`). The transport velocity used is the previous-iterate value \f$u^{\text{old}}\f$ (Picard linearisation). On non-orthogonal meshes the cell-tangential auxiliary DOFs enter via a face-to-cell projection that uses the second-order velocity reconstruction described below.

**Diffusive term** — \f$\nabla\cdot\bigl(\nu(\nabla u + (\nabla u)^\top)\bigr)\f$ is treated componentwise: each scalar component \f$u^d\f$ is fed to `fgrad` with the diffusivity \f$\nu\f$, and the resulting face flux feeds the viscous block. The reduction \f$\nabla\cdot\bigl(\nu(\nabla u + (\nabla u)^\top)\bigr) = \nu\Delta u\f$ holds only when \f$\nu\f$ is constant *and* \f$\nabla\cdot u = 0\f$ holds discretely; the implementation keeps the full Newtonian stress form so that variable viscosity is handled correctly.

### Cell-tangential velocity reconstruction

The viscous flux at a face needs the velocity vector at the adjacent cell barycenter — not only the face-normal component. The cell-tangential auxiliary DOFs \f$u_e^d\f$ carry this information, but they are **not** independent unknowns: after every linear solve, `Masse_PolyMAC_MPFA_Face::corriger_solution` overwrites them via a least-squares reconstruction from the face-normal values:

- **First order** (`Champ_Face_PolyMAC_MPFA::update_ve`, Whitney-1 pull-back): exact on constant fields; the cell-gradient it builds is \f$O(1)\f$ — not sufficient to make the viscous operator consistent.
- **Second order** (`Champ_Face_PolyMAC_MPFA::update_ve2`, default): exact on affine fields, obtained by solving a moment-fitting least-squares problem with LAPACK `DGELSY`. The runtime flag `INTERP_VE1` falls back to first order for debugging.

Second-order reconstruction is what makes the viscous operator consistent on smooth solutions. The soft mass-matrix entry on the cell-tangential rows is numerical-conditioning glue only — the *physical* closure is `update_ve2`.

### Pressure gradient

`Op_Grad_PolyMAC_MPFA_Face` reuses `fgrad` in pressure mode (`is_p = 1`). Two boundary regimes need care:

- On a **pressure-imposed** boundary (`Neumann_sortie_libre`), the trace value used in the multi-point formula is the prescribed value \f$p_b = p_{\text{imp}}(x_f)\f$. Accuracy of MPFA-O is preserved verbatim.
- On **walls and symmetry** boundaries, pressure has **no DOF and no externally-imposed value**. The wall-trace value \f$g_b(f) \approx p(x_f)\f$ is reconstructed from the discrete cell-tangential momentum equation of the wall-adjacent cell, by contracting it against \f$\vec{n}_f / (|f|\,|e|)\f$:

\f[
g_b(f) \;\approx\; \frac{1}{|f|\,|e|}\sum_d n_f(d)\,\Bigl(R^d_e \;-\; L^d_e(u) \;-\; \text{(interior-face terms)}\Bigr).
\f]

This `gb(f)` is split into a *constant part* (`gb`, from the explicit RHS) and a *velocity-dependent part* (`dgb_v_`) folded back into the velocity block of the matrix at the end of the assembly — so that the wall-pressure trace is treated **fully implicitly**, preserving the second-order accuracy of MPFA-O up to the wall. The naive alternative (homogeneous Neumann \f$\partial_n p = 0\f$) would degrade to first order on non-uniform meshes.

### Pressure-Poisson

Once `fgrad` is built, the discrete Poisson operator used by the projection step is

\f[
[\mathrm{div}\,\nabla p]_e \;=\; \sum_{f \in F_e} \sigma_{e,f}\,|f|\,\varphi_f\,[\vec{n}_f\cdot\nabla p]_f,
\f]

assembled in `Assembleur_P_PolyMAC_MPFA`. When the problem has no pressure reference (no `Neumann_sortie_libre`), the constant-pressure mode is pinned by doubling the (0, 0) entry on the first MPI rank carrying real cells.

## Boundary conditions

| BC type                        | Handled in                  | Effect on the local system                                |
|--------------------------------|-----------------------------|-----------------------------------------------------------|
| `Dirichlet`                    | `fgrad` (trace eq)          | hard-set \f$u_{f_v} = u_b\f$; RHS gets \f$\texttt{phif\_c}\cdot u_b\f$ |
| `Neumann`                      | `fgrad` (trace eq)          | RHS gets \f$+\|f_v\|\,\phi_{\text{imp}}\f$                |
| `Echange_impose`               | `fgrad` (trace eq)          | Robin-type coupling, RHS gets \f$h\,T_{\text{ext}}\f$     |
| `Frottement_*_impose`          | `fgrad` (trace eq)          | friction term on the cell or face-trace column            |
| `Symetrie` / outflow           | implicit                    | nothing to add (zero coefficient)                          |
| `Echange_contact_PolyMAC_MPFA` | `Op_Diff` helper            | vertex-level continuity across two problems via MEDCoupling; the only BC that produces off-problem matrix entries |

## PolyMAC_HFV

PolyMAC_HFV is based on a Hybrid Finite Volume (HFV) approach @cite EGH07 @cite DEGH10 and is mathematically close to PolyMAC_CDO, since HFV and CDO methods are equivalent @cite DEGH10. PolyMAC_MPFA differs from PolyMAC_HFV only in the face-flux reconstruction: the rest of the family (DOF layout, mass solver, convection operator, pressure assembly, multiphase support) is inherited from PolyMAC_HFV.
