@page Disc_VEF VEF

Initially introduced in @cite LM89, *Volume Élément Finis (VEF)* is a variant of the standard finite element and finite volume methods. The formalism developed in @cite Emonot1992 was subsequently used for the implementation of this method in **TRUST**.

## Finite Volume Element Method

### Core Idea

Consider the following instationary problem with velocity \f$\boldsymbol{u}\f$, flux term \f$\boldsymbol{F}\f$ and source term \f$\boldsymbol{S}\f$:

\f[
\partial_t \boldsymbol{u} + \nabla \cdot \boldsymbol{F} = \boldsymbol{S}
\f]

We introduce the control volume \f$\omega_f\f$ (see figure below) in which we evaluate the velocity \f$\boldsymbol{u}\f$, and integrate on \f$\omega_f\f$ between times \f$t^n\f$ and \f$t^{n+1}\f$:

\f[
\int_{\omega_f} (\boldsymbol{u}^{n+1} - \boldsymbol{u}^n)\mathrm{d}\boldsymbol{V}
+ \int_{\partial\omega_f} \int_{t^n}^{t^{n+1}} \boldsymbol{F} \cdot \boldsymbol{n} \, \mathrm{d}\boldsymbol{s}
= \int_{\omega_f} \int_{t^n}^{t^{n+1}} \boldsymbol{S} \, \mathrm{d}\boldsymbol{V}
\f]

For the Stokes equation: \f$\boldsymbol{F} = \mu \nabla \boldsymbol{u} - p\boldsymbol{I}\f$.
For the Navier-Stokes equation: \f$\boldsymbol{F} = \mu \nabla \boldsymbol{u} - p\boldsymbol{I} + \rho \boldsymbol{u} \otimes \boldsymbol{u}\f$.

![Control volume for velocity](figures/control_volume_velocity.png)

### Finite Volume Approach

Given a tetrahedral mesh \f$\mathcal{M}\f$, define points \f$\boldsymbol{x}_f\f$ as the barycenter of face \f$f\f$. The control volume \f$\omega_f\f$ links the vertices of face \f$f\f$ with the barycenters of the two cells \f$e_1, e_2\f$ sharing that face. Discretizing the evolution term:

\f[
\int_{\omega_f} \boldsymbol{u}^{m} \mathrm{d}\boldsymbol{V} \approx |\omega_f| \boldsymbol{u}_f^m \qquad m \in \{n, n+1\}
\f]

The flux term discretization yields:

\f[
|\omega_f|(\boldsymbol{u}_f^{n+1} - \boldsymbol{u}_f^n)
+ \Delta t^{n,n+1} |l_f| (\boldsymbol{F}^m_{e_2} - \boldsymbol{F}^m_{e_1})\,\vec{n}_f
= \Delta t^{n,n+1} \boldsymbol{S}_f^{n,n+1}
\f]

### Finite Element Basis

The VEF method uses the Crouzeix-Raviart basis: the full velocity vector is evaluated at the center of the faces of each cell, and within each cell the pressure is a constant at the cell center. Let \f$(\phi_f)_{f\in \mathcal{I}_F}\f$ be the velocity basis and \f$(\mathbb{I}_{e_k})_{k\in \mathcal{I}_E}\f$ the pressure basis:

\f[
\boldsymbol{u}_h = \sum_{f\in \mathcal{I}_F} \boldsymbol{u}_f \phi_f, \qquad
p_h = \sum_{k\in \mathcal{I}_E} p_k \mathbb{I}_{e_k}
\f]

![Control volumes for VEF-P0](figures/triangle.png)

### Discretization of the Flux Term (Stokes)

For the Stokes equation \f$\boldsymbol{F} = \mu \nabla \boldsymbol{u} - p\boldsymbol{I}\f$, the discrete gradient writes:

\f[
\int_{\partial\omega_f} \boldsymbol{\nabla} \phi_{f'} \cdot \boldsymbol{n} \, \mathrm{d}\boldsymbol{s}
= -\sum_{e \in \mathcal{M}} \frac{1}{|e|}\,\boldsymbol{S}_e^{f'} \cdot \boldsymbol{S}_e^f
\f]

and the pressure part:

\f[
\sum_{k \in \mathcal{I}_E} p_k \int_{\partial\omega_f \cap e_k} \boldsymbol{n} \, \mathrm{d}\boldsymbol{s}
= |l_f|(p_{e_2} - p_{e_1})\,\vec{n}_f
\f]

### Variational Formulation

Find \f$(\boldsymbol{u}_h, p_h) \in \mathbb{X}_h \times \overset{\circ}{\mathbb{N}}_h\f$ such that:

\f[
\left\{
\begin{aligned}
\partial_t m_h^V(\boldsymbol{u}_h,\boldsymbol{v}_h) + a_h^V(\boldsymbol{u}_h, \boldsymbol{v}_h) + b_h^V(\boldsymbol{v}_h, p_h) &= L_h^V(\boldsymbol{v}_h) & \forall \boldsymbol{v}_h \in \mathbb{X}_h, \\
c_h^V(\boldsymbol{u}_h, q_h) &= 0 & \forall q_h \in \overset{\circ}{\mathbb{N}}_h.
\end{aligned}
\right.
\f]

## Mathematical Properties

According to @cite Heib2003, the scheme satisfies the following properties:

- **Inf-sup condition**: ensures stability of the numerical scheme.
- **Continuity at edge midpoints**: implies weak continuity of velocity and enforces local mass conservation, leading to a divergence-free condition in each cell.
- **Well-posedness**: guarantees existence and uniqueness of the discrete solution.
- **Convergence rate for pressure**: order 1 in the \f$L^2\f$ norm.
- **Convergence rate for velocity**: order 2 in the \f$\boldsymbol{L^2}\f$ norm (convex domain).

Spurious currents for low velocities can appear when using the VEF approach @cite Fortin2006.

## Enriched Pressure Basis: VEF - \f$\mathbb{P}^{nc}/\mathbb{P}^0+\mathbb{P}^1\f$

To reduce spurious currents (useful for low viscosities), a pressure-enriched basis was studied in @cite Heib2003 @cite Fortin2006 and implemented in **TRUST** under the name VEF - \f$\mathbb{P}^{nc}/\mathbb{P}^0+\mathbb{P}^1\f$.

The idea is to add pressure unknowns \f$\mathbb{P}^1\f$ at the vertices of each cell, introducing a new control volume for the mass conservation. Two control volumes are used:

- \f$e_k\f$ for the constant part \f$\mathbb{P}^0\f$
- \f$\Pi_{v_i}\f$ for the \f$\mathbb{P}^1\f$ part associated with the vertex \f$v_i\f$

![Control volume for pressure P0 and P1](figures/pi_si_kl.png)

The stability of this basis is proved in @cite JCS23 and the inf-sup condition in @cite Fortin2006. This scheme is the most widely used VEF discretization in **TRUST**.
