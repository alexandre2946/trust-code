@mainpage TRUST Documentation

**TRUST** is a High Performance Computing (HPC) thermohydraulic engine for Computational Fluid Dynamics (CFD) developed at the Departement of %System and Structure Modelisation (DM2S) of the French Atomic Energy Commission (CEA).

The software was originally designed for conduction, incompressible single-phase, and Low Mach Number (LMN) flows with a robust Weakly-Compressible (WC) multi-species solver. A huge effort has been conducted recently, and now **TRUST** is able to simulate real compressible multi-phase flows.

**TRUST** is also progressively ported to support GPU acceleration, using the [Kokkos](https://kokkos.org/kokkos-core-wiki/) library. It has been selected to be a demonstrator of the [CExA](https://cexa-project.org/) project.

**TRUST** serves as the kernel of several CEA application codes. The [TrioCFD](https://triocfd.cea.fr/) software is an open source example. This software is OpenSource (**[BSD license](https://github.com/cea-trust-platform/trust-code/blob/master/License.txt)**), available **[here](https://github.com/cea-trust-platform/trust-code)** on GitHub.

---

- @subpage QuickStart
- @subpage UserGuide
- @subpage KeywordsReference
- @subpage UserTutorials
- @subpage TRUSTTools
- @subpage TRUST_Validation
- @subpage DevCorner
- @subpage References

---

<h2>Credits</h2>

**TRUST** relies on the following open source products: [SALOME](https://www.salome-platform.org/), [MEDCoupling](https://github.com/SalomePlatform/medcoupling), [ICoCo](https://github.com/cea-trust-platform/icoco-coupling), [VisIt](https://visit-dav.github.io/visit-website/index.html).

This page is powered by [Doxygen](https://www.doxygen.nl/) and [doxygen-awesome](https://github.com/jothepro/doxygen-awesome-css).
