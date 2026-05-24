@page Disc_PolyMAC_CDO PolyMAC_CDO

PolyMAC_CDO is based on **Compatible Discrete Operator (CDO)** schemes
@cite Bonelle2014 @cite Milani2020, designed to apply mimetic finite-volume
methods to non-conforming polyhedral meshes with the sole constraint that the
cells be star-shaped. In TRUST data files, the keyword `PolyMAC` is an alias
for `PolyMAC_CDO`.

The major distinction with PolyMAC_MPFA is that vorticity is an **explicit
unknown** of the velocity equation, not a derived quantity, and the viscous
term is discretised through a `curl–curl` decomposition (see "Diffusion
operator" below).

**Limitation.** PolyMAC_CDO does not support multiphase problems
(`Pb_Multiphase`). For multiphase, use `PolyMAC_MPFA` or `PolyMAC_HFV`.

## Dual mesh

PolyMAC introduces a dual mesh. The gravity center of each control volume
\f$cv \in \{E, F, \Sigma, V\}\f$ (cells, faces, edges, vertices) is called
\f$x_{cv}\f$. The dual entities are:

- **Dual cell** \f$\widetilde{v}\f$: located at the gravity center of the cell \f$x_e\f$.

![Dual cell when using PolyMAC_CDO](figures/Dual_cell.png)

- **Dual face** \f$\widetilde{\sigma}\f$: the line linking the gravity center of face \f$x_f\f$ to the gravity centers of the neighbouring cells.

![Dual face when using PolyMAC_CDO](figures/Dual_face.png)

- **Dual edge** \f$\widetilde{f}\f$: the surface linking the gravity centers of all neighbouring cells \f$x_e\f$, neighbouring faces \f$x_f\f$, and the edge gravity center \f$x_\sigma\f$.

![Dual edge when using PolyMAC_CDO](figures/Dual_edge.png)

## Location of the unknowns

| Field         | Location                                                                 |
|---------------|--------------------------------------------------------------------------|
| Pressure \f$p\f$   | Primal cell \f$e\f$ (= dual cell \f$\widetilde{v}\f$) **plus** an auxiliary face value \f$p_f\f$ |
| Velocity \f$u_f\f$ | Face-normal: \f$u_f = \frac{1}{\|f\|}\int \vec{u}\cdot\vec{n}_f\,\mathrm{d}S\f$         |
| Vorticity \f$\omega_\sigma\f$ | Edge-tangent in 3D: \f$\omega_\sigma = \frac{1}{\|\sigma\|}\int \vec{\omega}\cdot\vec{\tau}_\sigma\,\mathrm{d}l\f$; vertex scalar in 2D |

Scalar fields (e.g. temperature) carry a **cell DOF plus a face DOF**:
`Champ_Elem_PolyMAC_CDO` has \f$n_e + n_f\f$ values. The face DOF is the
"trace" unknown of the diffusion saddle-point (see below), exposed to
`Echange_contact_PolyMAC_CDO` for two-problem coupling. Pressure follows
the same layout.

## CDO scheme

The CDO scheme is based on exact operators that switch between localizations
on the primal and dual mesh.

![Paths between primal and dual mesh entities and corresponding operators](figures/CDO.png)

Operators are classified into two types:

- **Exact operators**: gradient, divergence, curl — linear operators encoded
  inline as signed incidence loops over the mesh connectivity, with no
  explicit matrix.
- **Approximated operators**: Hodge operators establish discrete closure
  relations between primal and dual meshes. They are represented by local
  matrices stored per cell.

### Exact discrete operators

**Discrete gradient:**

\f[
\text{GRAD}: V \rightarrow \Sigma, \qquad \forall \sigma \in \Sigma,\;
\text{GRAD}(\mathbf{p})|_{\tilde{\sigma}} = \sum_{\tilde{v} \in V_\sigma} \iota_{v,\sigma}\,p_v
\f]

**Discrete divergence:**

\f[
\text{DIV}: F \rightarrow E, \qquad \forall e \in E,\;
\text{DIV}(\phi)|_e = \sum_{f \in F_e} \iota_{f,e}\,\phi_f
\f]

**Discrete curl:**

\f[
\text{CURL}: \Sigma \rightarrow F, \qquad \forall f \in F,\;
\text{CURL}(\mathbf{u})|_f = \sum_{\sigma \in \Sigma_f} \iota_{\sigma,f}\,\mathbf{u}_\sigma
\f]

All operators have their dual counterparts \f$\widetilde{\text{GRAD}}\f$,
\f$\widetilde{\text{DIV}}\f$, \f$\widetilde{\text{CURL}}\f$ on the dual mesh.

### Hodge operators

Closure between primal and dual localisations uses Hodge operators
\f$H_{\alpha^{-1}}^{\widetilde{\mathcal{X}}\mathcal{Y}}\f$, expressed using
local geometry and material properties (diffusivity, viscosity). Local Hodge
operators must be **symmetric**, **locally stable** and
**\f$\mathbb{P}_0\f$-consistent** (exact on constant fields).

![Partitioning of the cell into sub-volumes attached to face (left) and edge (right)](figures/subvolume.png)

Three Hodge matrices are built once at domain-discretization time and reused
throughout the operators:

- \f$M_2\f$ (and its inverse \f$W_2\f$): face-normal velocity \f$\to\f$ tangent-to-dual-face, weighted by \f$\nu\f$. Built by `Domaine_PolyMAC_CDO::init_m2`.
- \f$M_1\f$: edge-tangent vorticity \f$\to\f$ integrated flux through the dual-face boundary. Built by `Domaine_PolyMAC_CDO::init_m1`.

For cell \f$e\f$ with \f$n_f\f$ faces, the local reconstruction is a dense
\f$n_f\times n_f\f$ matrix combining a **consistent (\f$\mathbb{P}_0\f$-exact)
part** with a **stabilisation correction** parameterised by

\f[
\beta =
\begin{cases}
1/D       & \text{on simplices (}n_f = D + 1\text{) — DGA reconstruction} \\
1/\sqrt{D} & \text{on general polyhedra — SUSHI reconstruction}
\end{cases}
\f]

The choice of \f$\beta\f$ is made **per cell** automatically by
`Domaine_PolyMAC_CDO::M2`.

#### Two paths for M2 / W2

`init_m2` dispatches between two algorithms via the global flag
`Option_PolyMAC_family::USE_NEW_M2`:

- **Default path** (`USE_NEW_M2 = 1`, `init_m2_new`): closed-form
  DGA/SUSHI formula. Stabilisation is baked into the formula via the
  \f$\beta\f$ correction; the build is deterministic and cannot fail.
- **Legacy path** (activated by `Option_PolyMAC use_osqp { }`, `init_m2_osqp`):
  solves a per-cell **quadratic-programming problem** via the OSQP library to
  find a W matrix that minimises off-diagonal \f$L_2\f$ norm subject to
  consistency \f$W \cdot R = N\f$ and bounded spectrum. Can produce sparser
  stencils but is more expensive and can fail on very distorted meshes
  (warning: `Matrice non stabilisee : attention le maillage est trop deforme!!`).

At first call the code prints a single diagnostic line summarising the
stabilisation result. The `% diag` / `% sym` / `lambda` fields are **only
populated under the OSQP path** — on the default path, only the `width`
(average stencil size per cell row) is informative.

## Scalar diffusion

For the elliptic equation
\f$-\underline{\nabla}\cdot(\underline{\underline{\kappa}}\,\underline{\nabla} T) = s\f$,
the CDO cell-based scheme yields the saddle-point system

\f[
\begin{pmatrix}
H_{\kappa^{-1}}^{\widetilde{F}\Sigma} & \widetilde{\mathbb{G}} \\
\mathbb{D} & 0
\end{pmatrix}
\begin{bmatrix} \phi \\ p \end{bmatrix}
= \begin{bmatrix} 0 \\ R_E(s) \end{bmatrix}
\f]

with \f$\widetilde{\mathbb{G}} = -\mathbb{D}^T\f$. The cell unknowns
\f$T_e\f$ satisfy the discrete divergence equation; the face trace unknowns
\f$T_f\f$ satisfy a continuity equation; the two blocks are coupled through
the \f$W_2\f$ Hodge matrix (`Op_Diff_PolyMAC_CDO_Elem`).

The trace value \f$T_f\f$ is **exact on affine fields** — this is the
\f$\mathbb{P}_0\f$-consistency property of the local Hodge.

A non-linear variant (`Op_Diff_Nonlinear_PolyMAC_CDO_Elem`) adds a
Le Potier / Mahamane-type monotonicity correction for high-Péclet regimes,
re-evaluated at each Newton iteration through cell- and face-side residual
indicators.

## Viscous diffusion of velocity (`curl–curl` form)

The viscous term in the momentum equation is **not** discretised directly as
\f$-\mu\Delta u\f$. Using the identity

\f[
-\underline{\Delta}\underline{u} \;=\; \underline{\nabla}\times\underline{\nabla}\times\underline{u} - \underline{\nabla}(\nabla\cdot\underline{u}) \;=\; \underline{\nabla}\times\underline{\omega}
\quad \text{(when }\nabla\cdot u = 0\text{)}
\f]

and introducing the vorticity \f$\underline{\omega} = \underline{\text{curl}}\,\underline{u}\f$
as an **explicit unknown** at the edges (3D) or vertices (2D), the cell-based
Stokes scheme @cite Bonelle2014 reads:

Find \f$(\mathbf{p}^*, \phi, \mathbf{\omega}^*) \in \widetilde{V} \times F \times \Sigma\f$ such that

\f[
\begin{pmatrix}
-H_{\mu^{-1}}^{\Sigma\widetilde{F}} & \widetilde{\mathbb{C}} \cdot H_{\rho^{-1}}^{F\widetilde{\Sigma}} & 0 \\
H_{\rho^{-1}}^{F\widetilde{\Sigma}} \cdot \mathbb{C} & 0 & \widetilde{\mathbb{G}} \\
0 & -\mathbb{D} & 0
\end{pmatrix}
\begin{bmatrix} \mathbf{\omega}^* \\ \phi \\ \mathbf{p}^* \end{bmatrix}
= \begin{bmatrix} 0 \\ S(\rho, \underline{f}^*) \\ 0 \end{bmatrix}
\f]

The edge rows enforce the **definition of vorticity**
(\f$M_1 \cdot \omega = \nu\,\text{CURL}(v)\f$); the face rows enforce the
**viscous flux balance** with \f$M_2 \cdot \text{CURL}(\omega)\f$ on the
left-hand side. The unknowns of the velocity field are
\f$(u_f)_{f \in F} \cup (\omega_\sigma)_{\sigma \in \Sigma}\f$ — total
\f$n_f + n_\sigma\f$ DOFs.

The mass matrix on velocity is the \f$M_2\f$ Hodge itself — it is **not
diagonal**. On steady-state problems with fixed mesh, the factorisation is
built once via PETSc/Cholesky and reused
(`Domaine_PolyMAC_CDO::init_m2solv`).

## Pressure gradient and divergence

The pressure gradient on faces is the simple two-point difference

\f[
[\vec{\nabla} p]_f \cdot \vec{n}_f \, |f| \;=\; |f|\,\varphi_f\,(p_{e_2} - p_{e_1})
\f]

(`Op_Grad_PolyMAC_CDO_Face`), and the divergence on cells is

\f[
[\nabla\cdot u]_e \;=\; \sum_{f \in F_e} \sigma_{e,f}\,|f|\,\varphi_f\,u_f
\f]

(`Op_Div_PolyMAC_CDO`). The pair
\f$(\mathbb{D}, \widetilde{\mathbb{G}} = -\mathbb{D}^\top)\f$ discretises the
saddle-point block of the incompressible NS system. The natural primal-dual
orthogonality of CDO replaces the multi-point reconstruction needed by
PolyMAC_MPFA; in particular, **there is no wall-pressure trace closure** —
walls and symmetry faces do not contribute to `Op_Grad`, and the right
boundary closure is enforced through the \f$M_2\f$ Hodge and the curl-curl
viscous term.

## Pressure-Poisson assembly

The pressure-Poisson matrix used by the projection step
(`Assembleur_P_PolyMAC_CDO`) has size \f$(n_e + n_f)^2\f$: **face pressures
participate as auxiliary unknowns**, derived from cell pressures through the
\f$W_2\f$ Hodge during the pressure-solver step. For each cell \f$e\f$, the
local contribution is

\f[
M_{f f_b} \;=\; \frac{|f|\,|f_b|\,\varphi_e}{|e|}\,W_2(i, j)
\f]

assembled on the face-face block, with the corresponding \f$-M_{e f}\f$
entries on the cell-face cross-blocks.

- A `Neumann_sortie_libre` BC fixes the face-pressure equation to
  \f$\delta p_f = 0\f$ and reinjects the imposed value through the RHS.
- In the absence of any pressure-Dirichlet BC, the constant-pressure mode is
  pinned by doubling the (0, 0) entry on the first MPI rank carrying real cells.

## Convection

### Scalar convection

`Op_Conv_Amont_PolyMAC_CDO_Elem` and `Op_Conv_Centre_PolyMAC_CDO_Elem` are
the classical first-order upwind and centred schemes on face-normal
velocities:

\f[
\text{flux}_f =
\begin{cases}
\Delta t \cdot u_f \cdot |f|\,\varphi_f \cdot T_{e_{\text{up}}} & \text{(upwind)} \\
\Delta t \cdot u_f \cdot |f|\,\varphi_f \cdot \frac{T_{e_1} + T_{e_2}}{2} & \text{(centred)}
\end{cases}
\f]

The auxiliary face trace \f$T_f\f$ introduced by the diffusion operator is
**not** used in convection — only cell values.

### Velocity convection

`Op_Conv_EF_Stab_PolyMAC_CDO_Face` is a face-based blended EF-Stab scheme:

\f[
\text{dfac}(j) \;=\; \frac{|f|\,u_f\,\varphi_e}{2}\,\bigl(1 + \alpha\,\mathrm{sgn}(u_f)\bigr),\qquad \alpha \in [0, 1]
\f]

with \f$\alpha\f$ blending pure upwind (\f$\alpha = 1\f$, exposed as
`Op_Conv_Amont_PolyMAC_CDO_Face`) and pure centred (\f$\alpha = 0\f$,
`Op_Conv_Centre_PolyMAC_CDO_Face`). Two coupling paths:

- **Face-equivalence path** (typical on K-orthogonal meshes): when an
  equivalent face exists across the cell, contributions are assembled
  directly face-to-face with a porosity / normal-sign multiplier. This is the
  path that recovers MAC on Cartesian meshes.
- **Element-decomposition path** (general polyhedra): cousin faces of the
  two cells touching \f$f\f$ are coupled via the geometric inner product
  \f$(\vec{x}_{f_c} - \vec{x}_{e_b}) \cdot (\vec{x}_{f_b} - \vec{x}_e)\f$,
  projected through the face normals.

This operator is **single-phase only** (`assert(N == 1)` in the code).

There is **no F → cell reconstruction step inside the convection operator**.
The Perot face-to-cell formula

\f[
\vec{u}_e \;=\; \frac{1}{|e|}\sum_{f \in F_e} \sigma_{e,f}\,|f|\,(\vec{x}_f - \vec{x}_e)\,u_f
\f]

@cite Bonelle2014 (and Perot–Basumatary, J. Comput. Phys. 2014) **is**
implemented in `Domaine_PolyMAC_CDO::init_ve`, but it is used only for
output / post-processing, pressure-loss source terms, and the boundary CURL
handling in the viscous operator — not inside the convective stencil.

## Boundary conditions

For scalars (`Champ_Elem_PolyMAC_CDO`):

| BC type                          | Effect on the local system                                       |
|----------------------------------|------------------------------------------------------------------|
| `Dirichlet`                      | hard-set \f$T_f = T_{\text{imp}}\f$; identity row on the face DOF |
| `Neumann_paroi`                  | RHS gets \f$+\|f\|\,\phi_{\text{imp}}\f$                          |
| `Echange_externe_impose`         | RHS gets \f$h\,T_{\text{ext}}\f$; matrix gets \f$+\|f\|\,h\f$ on (f, f) |
| `Echange_global_impose`          | adds friction-type term on the cell column                       |
| `Echange_contact_PolyMAC_CDO`    | monolithic two-problem coupling via \f$W_2\f$ stencil exchange — the only BC that produces matrix entries across two problems |
| `Symetrie` / outflow             | nothing to add                                                   |

For velocity (`Champ_Face_PolyMAC_CDO`):

| BC type                                   | Effect on the local system                                          |
|-------------------------------------------|---------------------------------------------------------------------|
| `Dirichlet` (\f$\vec{u}\f$ imposed)       | \f$u_f\f$ hard-set in mass row; full \f$\vec{u}\f$ vector fed to the boundary CURL via `ra` |
| `Symetrie` (\f$\vec{n}_f\cdot\vec{u} = 0\f$) | \f$u_f = 0\f$; tangential \f$\omega\f$ handled through the `va` reconstruction |
| `Neumann_sortie_libre` (pressure imposed) | pressure Dirichlet; \f$u_f\f$ free                                  |

## Incompressible Navier-Stokes

Assembling the building blocks above gives, at each Newton step, the discrete
incompressible NS system

\f[
\begin{cases}
H_{\mu^{-1}}\,\underline{\omega}^* = \underline{\nabla}\times(\rho^{-1}\underline{\phi}), \\
\partial_t(\rho^{-1}\underline{\phi}) + \nabla\cdot((\rho^{-1}\underline{\phi})\otimes(\rho^{-1}\underline{\phi}))
+ \rho^{-1}\underline{\nabla}\times\underline{\omega}^* + \underline{\nabla}p^* = \underline{f}^*, \\
\underline{\nabla}\cdot\underline{\phi} = 0.
\end{cases}
\f]

The continuous-to-discrete correspondence is summarised by:

| Continuous identity                                              | Discrete realisation                                  |
|------------------------------------------------------------------|-------------------------------------------------------|
| \f$-\Delta u = \nabla\times\omega\f$ (when \f$\nabla\cdot u = 0\f$) | \f$M_2 \cdot \text{CURL}(\omega) = \nu \cdot M_1^{-1} \cdot \text{CURL}(u)\f$ |
| \f$\omega = \nabla\times u\f$                                    | \f$M_1\,\omega = \nu\,\text{CURL}(u)\f$ (edge row)    |
| \f$\nabla\cdot u = 0\f$                                          | \f$D\,u_f = 0\f$ (cell row)                           |
| \f$\nabla p\f$                                                   | \f$-D^\top p_e\f$ (face row contribution)             |
| \f$\partial_t u\f$                                               | \f$M_2 / \Delta t \cdot u_f\f$ (mass relation)        |

This is the **discrete Hodge–Helmholtz decomposition** built into the
geometry of the dual mesh: the CDO scheme preserves it exactly, with
\f$D \cdot \widetilde{G} = D \cdot (-D^\top) \equiv\f$ discrete pressure
Laplacian.

## Practical diagnostics

- The first useful diagnostic at startup is the `init_m2` line. Under the
  default `USE_NEW_M2 = 1`, only the `width` field is meaningful (average
  stencil size); under `Option_PolyMAC use_osqp { }`, the `% diag`, `% sym`
  and eigenvalue-range fields are populated too.
- `Matrice non stabilisee : attention le maillage est trop deforme!!`
  (legacy OSQP path only) signals that the QP could not produce a positive-definite
  \f$W\f$ — convergence is at risk; consider remeshing or switching back to
  the default path.

## PolyMAC_HFV

PolyMAC_HFV is based on a Hybrid Finite Volume approach @cite EGH07 @cite DEGH10.
It is mathematically close to PolyMAC_CDO since HFV and CDO methods are
equivalent @cite DEGH10, but implementation-wise PolyMAC_HFV shares the
architecture of PolyMAC_MPFA (only the face-flux reconstruction differs).
