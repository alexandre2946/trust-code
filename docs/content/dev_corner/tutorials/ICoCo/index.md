@page DevTut_Icoco TRUST ICoCo Tutorial

This tutorial introduces the ICoCo coupling interface as implemented in the TRUST
platform. It covers the full workflow from environment setup to bidirectional field
exchange and Python-based coupling supervision.

The tutorial is structured as three progressive exercises:

- @subpage DevTut_Icoco_basics — Environment setup and first steps with a single-problem test case.
- @subpage DevTut_Icoco_coupled — Setting up and running a genuinely coupled problem with field exchange.
- @subpage DevTut_Icoco_python — Using the ICoCo API from a Python supervisor script.

---

## Why Code Coupling?

Traditional numerical simulation codes are designed around a single physical domain:
one code handles thermics, another handles mechanics, and so on. Real-world
engineering studies, however, require the simultaneous simulation of several coupled
physical phenomena. Nuclear reactor simulations, for example, must account for
thermics, neutronics, and structural mechanics in a mutually consistent fashion.

The natural answer to this requirement is **code coupling**: allowing several
specialized codes to communicate with one another while each remains responsible for
its own physics domain.

## The Role of the Supervisor

Any coupled computation requires a central entity to drive the overall calculation:
the **supervisor**. The supervisor is responsible for initializing all participating
codes, advancing them through their respective time loops, and centralizing the data
exchanges and unit conversions between them. In practice, the supervisor is typically
implemented as a C++ program or a Python script.

In a naive approach, the supervisor must be aware of the internal API of every code
it orchestrates. This quickly becomes unmanageable as soon as more than two codes are
involved, or when a single supervisor is intended to operate with different pairs or
groups of codes.

## A Common Interface

The standard solution is to define a **unique coupling interface** to which every
participating code must conform. The interface acts as a contract: it specifies a
fixed set of methods with standardized signatures, while leaving the implementation
details entirely to each individual code. The supervisor then only needs to know this
common interface, regardless of which specific codes are being coupled.

This design provides **versatility**: any code that correctly implements the interface
can be substituted into an existing supervisor with minimal or no modification. The
trade-off is that the interface must be implemented once for each new code that is to
be coupled.

## ICoCo: Interface for Code Coupling

**ICoCo** (Interface for Code Coupling) is the standard coupling interface used
within the TRUST ecosystem. Its key characteristics are the following.

It is written in C++ and, since TRUST 1.8.3, also wrapped via SWIG to allow use
from Python scripts. It is designed specifically for simulation codes structured
around iterative time loops. It defines a set of standard methods with fixed
signatures and no default implementation; notable methods include
`initializeTimeStep()`, `validateTimeStep()`, and `abortTimeStep()`. Data exchange
between codes is performed through the notion of a **field** — a set of scalar or
vector values defined on a mesh. The standard implementation relies on
`MEDCouplingFieldDouble` from the SALOME MEDCoupling library.

@note The ICoCo interface is formally described in the report
*"An Interface for Code Coupling, ICoCo v1.2"*, accessible within any TRUST
installation at `$TRUST_ROOT/doc/Kernel/ICoCo_V1.2.pdf`. The detailed API
specification ("APIProblem.pdf") can be requested from the TRUST support team
at trust@cea.fr.