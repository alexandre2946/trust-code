@page NM_BoundaryConditions Boundary conditions

@note This page is currently under construction. Additional boundary conditions will be added in future updates.

## Robin Boundary Conditions

@note Robin boundary conditions are only available with the **VEF Pnc/P0** discretization.

Robin boundary conditions consist of a linear combination between the flux term \f$\boldsymbol{F}\f$ and the variable term. They are useful for fluid-structure interactions or domain decomposition. The implementation in **TRUST** decomposes the condition into normal and tangential parts.

Let \f$\boldsymbol{n} = (n_x, n_y)\f$ be the outward normal vector. The flux term is:

\f[
\boldsymbol{F} = F_n\boldsymbol{n} + \boldsymbol{F_t}
\f]

For the Navier-Stokes equations:

\f[
\begin{aligned}
F_n &= \nu\nabla_n\boldsymbol{u}\cdot\boldsymbol{n} + \chi(\boldsymbol{u}\cdot\boldsymbol{n})(\boldsymbol{u}\cdot\boldsymbol{n}) - p \\
\boldsymbol{F_t} &= \nu\nabla_n\boldsymbol{u}\times\boldsymbol{n} + \chi(\boldsymbol{u}\cdot\boldsymbol{n})(\boldsymbol{u}\times\boldsymbol{n})
\end{aligned}
\f]

with \f$\nu\f$ the viscosity and \f$\chi \in \{0,1\}\f$. In 2D, the cross product is replaced by a projection onto the tangential vector \f$\boldsymbol{t} = (-n_y, n_x)\f$.

Two Robin parameters are defined: \f$\alpha\f$ for the normal part and \f$\beta\f$ for the tangential part, with Robin data:
- a normal scalar function \f$g_N\f$
- a tangential function \f$\boldsymbol{g_T}\f$ (scalar in 2D, vector in 3D)

The Robin boundary conditions implemented in **TRUST** are:

\f[
\begin{aligned}
\alpha F_n + \boldsymbol{u}\cdot\boldsymbol{n} &= g_N \\
\beta\boldsymbol{F_t} + \boldsymbol{u}\times\boldsymbol{n} &= \boldsymbol{g_T}
\end{aligned}
\f]

### Keyword: `Robin_VEF`

Parameters:
- `alpha`, `beta` as defined above.
- `champ_front_normal_et_tangentiel` followed by the field data (concatenation of \f$g_N\f$ and \f$\boldsymbol{n}\times\boldsymbol{g_T}\f$).

**2D example** — for \f$\boldsymbol{u}=(y,-x)\f$, \f$p=0.5(x^2+y^2)-1/3\f$, \f$\boldsymbol{n}=(1,0)\f$:

```
Robin_VEF {
    alpha 3
    beta 4
    champ_front_normal_et_tangentiel_robin champ_front_fonc_txyz 2
        -1.5*x^2-4.5*y^2+y+1.0
        4*x*y-x-4
}
```

The first function is \f$g_N\f$ and the second is \f$g_T\f$.

**3D**: the field `champ_front_normal_et_tangentiel` will have 4 components (one for \f$g_N\f$ and three for \f$\boldsymbol{n}\times\boldsymbol{g_T}\f$).

@note The notation \f$\boldsymbol{n}\times\boldsymbol{g_T}\f$ is used because we want to write the real tangential component of \f$\boldsymbol{u}\f$.
