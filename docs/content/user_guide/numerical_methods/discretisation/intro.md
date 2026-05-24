@page Disc_Intro Introduction

In the **TRUST** code, different numerical schemes are available to the user: VDF, VEF and the PolyMAC family.

- The VDF discretization is based on the Marker and Cell scheme @cite HW65.

- The VEF discretization is based on the Crouzeix-Raviart element method @cite Emonot1992.

The PolyMAC discretization family has been developed since 2018. Three PolyMAC variants are available in **TRUST**. They have been built using a Finite Volume (FV) framework on a staggered mesh so as to extend the MAC scheme @cite HW65 to complex grids:

- **PolyMAC**: based on a Compatible Discrete Operator (CDO) approach @cite Bonelle2014 @cite Milani2020.

- **PolyMAC_MPFA**: based on a Multi Point Flux Approximation (MPFA) approach @cite AM08 @cite D14 @cite LePotier2017.

- **PolyMAC_HFV**: based on a Hybrid Finite Volume (HFV) approach @cite EGH07 @cite DEGH10.

For each method the core ideas and the main steps for the discretization of the incompressible Navier-Stokes equations are presented. For now, the PolyMAC and PolyMAC_MPFA parts are completed; the others are a work in progress.

## Notations

Let's consider a space \f$\Omega\f$ and a certain grid \f$\mathcal{M}\f$ of non-overlapping polyhedra that map \f$\Omega\f$.

In the following:

- A polyhedron of the grid will be called a **cell**: \f$e\f$.
- A **face** \f$f\f$ is defined as the intersection of two cells or one cell and a boundary. Faces are supposed to be planar.
- An **edge** \f$\sigma\f$ is defined as the intersection of faces or faces and a boundary (3D only).
- A **vertex** \f$v\f$ is defined as the intersection of edges or edges and a boundary.

The set of cells is called \f$E\f$. The set of faces of cell \f$e\f$ is denoted \f$F_e\f$. The set of edges of face \f$f\f$ is \f$\Sigma_f\f$ and the two vertices of edge \f$\sigma\f$ are \f$V_{\sigma}\f$.

When a face \f$f\f$ separates two cells they are denoted \f$e_1\f$ and \f$e_2\f$ (the choice between the two is fixed by the orientation of \f$f\f$). The gravity centers (barycenters) of these entities are written \f$x_e\f$, \f$x_f\f$, \f$x_\sigma\f$, \f$x_v\f$.

The measure of an unknown \f$x\f$ at a control volume \f$cv\f$ is:

\f[
[x]_{cv} = \frac{1}{|cv|} \int_{cv} x \, \mathrm{d}(cv)
\f]

where \f$|\cdot|\f$ is a global measure operator: \f$|e|\f$ refers to the volume of cell \f$e\f$, \f$|f|\f$ to the surface of face \f$f\f$, and \f$|\sigma|\f$ to the length of edge \f$\sigma\f$.

For face-based fluxes we use the outward unit normal \f$\vec{n}_f\f$ together with the orientation sign

\f[
\sigma_{e,f} \;=\; \begin{cases} +1 & \text{if } \vec{n}_f \text{ points outward of cell } e, \\ -1 & \text{otherwise.} \end{cases}
\f]

(The symbol \f$\sigma\f$ is overloaded with the edge index — this is consistent with both the FV and CDO literature, and the context always makes the meaning clear.)

The main unknowns are \f$u\f$ (velocity) and \f$p\f$ (pressure); for scalar transport we use \f$T\f$ (temperature) and \f$s\f$ for source terms. When porosities are used we write \f$\varphi_e\f$ for the cell porosity and \f$\varphi_f\f$ for the face porosity.

**Vector typography.** Scalars are typeset plain (\f$u, p, T\f$). Vectors are decorated with an arrow when the distinction matters (\f$\vec{n}_f, \vec{x}_f, \vec{u}\f$). The VEF page uses bold (\f$\boldsymbol{u}, \boldsymbol{F}\f$) instead of arrows for legibility in operator-heavy formulas — this is the only place where the convention differs.
