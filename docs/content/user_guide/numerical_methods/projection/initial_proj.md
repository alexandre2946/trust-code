@page Proj_InitialProjection Description of the initial projection

In TRUST, an initial pressure is not required from the user. The determination of the initial pressure is handled in `Navier_Stokes_std::preparer_calcul`, using an approach similar to the correction step of the Chorin algorithm.

## Goal

We need to find the initial pressure that corrects the initial velocity to be divergence-free at every cell. Even if the initial velocity is globally divergence-free, it is not necessarily divergence-free at each cell of the mesh.

Let \f$U^{ini}\f$ be the initial velocity given by the user and \f$U^0\f$ the real velocity used to launch the algorithm. The velocity \f$U^0\f$ satisfies:

\f[
\begin{aligned}
\mathbb{M}\frac{U^0 - U^{ini}}{\delta t^0} + \mathbb{B}^t P^0 &= 0 \\
\mathbb{B} U^0 &= 0
\end{aligned}
\f]

## Initial Pressure Computation

Multiplying the first equation by \f$\mathbb{B}\mathbb{M}^{-1}\f$ and using the divergence-free condition, the initial pressure is computed by solving:

\f[
\delta t^0\,\mathbb{B}\mathbb{M}^{-1}\mathbb{B}^t P^0 = \mathbb{B} U^{ini}
\f]

The initial velocity is then recovered as:

\f[
U^0 = U^{ini} - \delta t^0\,\mathbb{M}^{-1}\mathbb{B}^t P^0
\f]

@note Some keywords allow modification of this initial pressure system: it is possible to impose an initial pressure or use the source term to find the correct pressure.

Once the initial values of velocity and pressure are determined, the projection algorithms can proceed. See @ref Proj_AvailableMethods for the list of available projection methods.
