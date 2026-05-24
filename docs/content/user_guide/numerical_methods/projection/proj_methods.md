@page Proj_AvailableMethods List of the available projection methods

## Chorin Algorithm

The first projection method was proposed by @cite C67 and @cite T69. The pressure is not included in the intermediate velocity prediction.

### Step 1: Velocity Prediction

Find intermediate velocity \f$U^*\f$ satisfying:

\f[
\mathbb{M}\frac{U^* - U^n}{\delta t^n} + \mathbb{A}U^* = F^{n+1}
\f]

The boundary conditions of \f$U^*\f$ are the same as those of \f$U^{n+1}\f$. Note that \f$U^*\f$ has no reason to be divergence-free.

### Step 2: Pressure Correction

Find \f$P^{n+1}\f$ from the following system (convection and Laplacian matrices are omitted as the mass matrix is easy to invert):

\f[
\begin{aligned}
\mathbb{M}\frac{U^{n+1} - U^*}{\delta t^n} + \mathbb{B}^T P^{n+1} &= 0 \\
\mathbb{B} U^{n+1} &= 0
\end{aligned}
\f]

Left-multiplying by \f$\mathbb{B}\mathbb{M}^{-1}\f$ yields a Poisson system on pressure:

\f[
\delta t^n\,\mathbb{B}\mathbb{M}^{-1}\mathbb{B}^T P^{n+1} = \mathbb{B}U^*
\f]

### Step 3: Update

\f[
U^{n+1} = U^* + \delta t^n\,\mathbb{M}^{-1}\mathbb{B}^T P^{n+1}
\f]

Summing steps 1 and 2, the velocity \f$U^*\f$ appears in the reconstructed system behind the convection and Laplacian — introducing a *splitting error* @cite G96.

@note In TRUST, this method is used for explicit schemes.

---

## Chorin Algorithm with Pressure Increment

A modification proposed by @cite G79 takes the pressure at the previous time in step 1 and solves for a pressure increment in step 2. Used for semi-implicit or implicit schemes.

### Step 1: Velocity Prediction

\f[
\mathbb{M}\frac{U^* - U^n}{\delta t^n} + \mathbb{A}U^* + \mathbb{B}^T P^n = F^{n+1}
\f]

### Step 2: Pressure Correction

Solve for the pressure increment \f$\delta P := P^{n+1} - P^n\f$:

\f[
\delta t^n\,\mathbb{B}\mathbb{M}^{-1}\mathbb{B}^T \delta P = \mathbb{B}U^*
\f]

### Step 3: Update

\f[
\begin{aligned}
U^{n+1} &= U^* + \delta t^n\,\mathbb{M}^{-1}\mathbb{B}^T \delta P \\
P^{n+1} &= P^n + \delta P
\end{aligned}
\f]

---

## SIMPLE Algorithm

Semi-Implicit Method for Pressure Linked Equations @cite PS72 @cite CGPS07. Similar to Chorin with pressure increment, but the mass matrix \f$\mathbb{M}/\delta t^n\f$ in steps 2 and 3 is replaced by the diagonal of the convection-diffusion matrix:

\f[
\mathbb{D} := \text{diag}\!\left(\mathbb{A} + \frac{\mathbb{M}}{\delta t^n}\right)
\f]

**Pressure correction:**

\f[
\mathbb{B}\mathbb{D}^{-1}\mathbb{B}^T \delta P = \mathbb{B}U^*
\f]

**Update:**

\f[
\begin{aligned}
U^{n+1} &= U^* + \mathbb{D}^{-1}\mathbb{B}^T \delta P \\
P^{n+1} &= P^n + \delta P
\end{aligned}
\f]

A relaxation can be applied to the pressure (or velocity) at the update step. Implementation details: `Simple.h`.

---

## SIMPLER Algorithm

SIMPLE Revised @cite P80 — adds a pre-computed pressure step before the prediction step.

### Step 0: Pre-compute the Pressure

Define \f$\mathbb{E} := \mathbb{A} + \mathbb{M}/\delta t^n - \mathbb{D}\f$. Find intermediate velocity \f$U^p\f$:

\f[
\mathbb{D}(U^n)\,U^p - \mathbb{E}U^n = F^{n+1}
\f]

Then solve for the pre-computed pressure:

\f[
\mathbb{B}\mathbb{D}^{-1}\mathbb{B}^t P^{n+1} = \mathbb{B}U^p
\f]

### Step 1: SIMPLE on \f$(U^{n+1}, P^{n+1})\f$

Apply the standard SIMPLE algorithm (velocity prediction → pressure correction → update). Implementation details: `Simpler.h`.

---

## PISO Algorithm

Pressure-Implicit with Splitting of Operators @cite I83 — a two-step projection method extending SIMPLE with a second correction that considers the non-diagonal part of \f$\mathbb{A}\f$.

### Step 1: SIMPLE Step

**Velocity prediction:**

\f[
\mathbb{M}\frac{U^* - U^n}{\delta t^n} + \mathbb{A}U^* + \mathbb{B}^T P^n = F^{n+1}
\f]

**First pressure increment:**

\f[
\mathbb{B}\mathbb{D}^{-1}\mathbb{B}^T \delta P^{p1} = \mathbb{B}U^*
\f]

**First update:**

\f[
\begin{aligned}
U^{p1} &= U^* + \mathbb{D}^{-1}\mathbb{B}^T \delta P^{p1} \\
P^{p1} &= P^n + \delta P^{p1}
\end{aligned}
\f]

### Step 2: Second Pressure Correction

The second correction considers the non-diagonal part \f$\mathbb{E}_A\f$ of the convection-diffusion matrix:

\f[
\mathbb{B}\mathbb{D}^{-1}\mathbb{B}^t \delta P^{p2} = \mathbb{B}\mathbb{D}^{-1}\mathbb{E}_A U^{p1}
\f]

**Final update:**

\f[
\begin{aligned}
U^{n+1} &= \mathbb{E}_A U^{p1} - \mathbb{B}^t \delta P^{p2} \\
P^{n+1} &= P^{p1} + \delta P^{p2}
\end{aligned}
\f]

Algebraic details are presented in @cite I83 or in `Piso.h`.
