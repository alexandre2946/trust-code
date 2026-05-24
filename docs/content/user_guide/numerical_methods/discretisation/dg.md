@page Disc_DG Discontinuous Galerkin Methods

Discontinuous Galerkin (DG) methods form a class of finite element methods particularly suited for solving partial differential equations. Unlike continuous Galerkin methods, DG allows for discontinuities between elements, providing greater flexibility and robustness for complex geometries or highly dynamic phenomena @cite HW07 @cite EDP12 @cite CKS12.

| **Pros** | **Cons** |
|---|---|
| High-order accuracy | Large number of unknowns |
| Handles non-conforming meshes | Numerous parameters to tune |
| High arithmetic intensity | |

A Symmetric Interior Penalty (SIP) method has been implemented in **TRUST** for solving the Non-Stationary Heat Equation. The SIP method is more performant than mixed DG methods for first- and second-order approximations.

## Definitions

Considering a face \f$f\f$ shared by two cells \f$e_1\f$ and \f$e_2\f$, the **interface average** of a quantity \f$y\f$:

\f[
\{\{ y \}\}_f(x) = \frac{1}{2}\left(y|_{e_1}(x) + y|_{e_2}\right)
\f]

The **interface jump** (when the normal of \f$f\f$ is defined from \f$e_1\f$ to \f$e_2\f$):

\f[
{\left[\!\left[ y \right]\!\right]}_f(x) = y|_{e_1}(x) - y|_{e_2}
\f]

and otherwise:

\f[
{\left[\!\left[ y \right]\!\right]}_f(x) = y|_{e_2}(x) - y|_{e_1}
\f]

![Definition of the average and jump notations](figures/scheme_jump_average.png)

## SIP DG Method for the Poisson Problem

Find \f$u \in H^1_0(\Omega)\f$ such that:

\f[
-\text{div}(k\nabla u) = s \quad \Rightarrow \quad a_{dg}(u_h, v_h) = \int_\Omega s\,v_h, \quad \forall v_h \in X_{DG}
\f]

**Discrete bilinear form:**

\f[
\begin{aligned}
a_{dg}(u_h, v_h) &= \sum_{e \in E} k_e \int_e \nabla u_h \cdot \nabla v_h \\
&\quad - \sum_{f\in F_e} \int_f k_f \{\{\nabla u_h\}\}_f \cdot \vec{n}_f\,\left[\!\left[v_h\right]\!\right] \\
&\quad - \sum_{f\in F_e} \int_f k_f \{\{\nabla v_h\}\}_f \cdot \vec{n}_f\,\left[\!\left[u_h\right]\!\right] \\
&\quad + \sum_{f\in F_e} \frac{\eta}{h_e} \int_f \left[\!\left[u_h\right]\!\right]_f\,\left[\!\left[v_h\right]\!\right]_f
\end{aligned}
\f]

where \f$h_e\f$ is the diameter of the circumscribed circle of \f$e\f$.

- First term: **consistency**
- Second and third terms: **symmetry**
- Last term: **stability**

The global stiffness matrix \f$\mathbf{K}\f$ has a block-structured form reflecting the element-wise stencil:

\f[
\mathbf{K} = \begin{bmatrix}
\mathbf{K}_{1,1} & \mathbf{K}_{1,2} & 0 & \cdots \\
\mathbf{K}_{1,2}^e & \mathbf{K}_{2,2} & \mathbf{K}_{2,3} & \cdots \\
0 & \mathbf{K}_{2,3}^e & \mathbf{K}_{3,3} & \cdots \\
\vdots & & & \ddots
\end{bmatrix}
\f]

![Possible mesh with the Discontinuous Galerkin discretization](figures/mesh_DG.png)

The stability parameter \f$\eta\f$ is computed automatically to ensure coercivity.

## Non-Stationary Heat Equation

For all \f$t \in [0, t_{max}]\f$, find \f$T(t) \in H^1_0(\Omega)\f$ such that:

\f[
\rho C_p \frac{dT}{dt} - \text{div}(k\nabla T) = s
\f]

| Quantity | Description |
|---|---|
| \f$k\f$ | Thermal conductivity (W·m⁻¹·K⁻¹) |
| \f$\rho\f$ | Density (kg·m⁻³) |
| \f$C_p\f$ | Heat capacity (J·kg⁻¹·K⁻¹) |
| \f$T\f$ | Temperature (K) |
| \f$s\f$ | Heat source (W·m⁻³) |

**Weak form:**

\f[
\rho C_p\,m_{dg}\!\left(\frac{dT_h}{dt},\theta_h\right) + a_{dg}(T_h,\theta_h)
= \int_\Omega s\,\theta_h, \quad \forall\theta_h \in X_{DG},\; \forall t\in[0,t_{max}]
\f]

**Time integration:**

- **Implicit Euler**: allows larger time steps; requires solving a linear system at each step.
- **Explicit Euler**: fast iterations; requires small time steps (stability constraint).

## DG Options

In your data file, add an `Option_DG` block:

```
Option_DG
{
    order 2
    gram_schmidt 1
}
```

`order` sets the discretization order (only **1 and 2** are currently available). `gram_schmidt 1` enables Gram-Schmidt orthonormalization of basis functions, which diagonalizes the mass matrix — useful with explicit schemes.
