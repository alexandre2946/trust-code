@page Disc_VDF VDF

The Finite Volume Difference (VDF) discretization (class `VDF_discretisation`, alias `VDF`) is the simplest and most efficient discretization of the TRUST platform. It is compatible with conforming meshes with hexahedral elements. Note: do not confuse a hexahedral mesh with a Cartesian mesh — the VDF mesh is not structured and does not follow the IJK indexing.

As stated by its name, the VDF is a conservative finite volume scheme of Marker-and-Cell (MAC) type @cite HW65. The discretization of each term of the equation is performed by integrating over a control volume. Diffusion gradient terms are approximated by a linear difference equation. All scalars are stored at the center of each control volume except the velocity field, which is defined on a staggered mesh.

This discretization **supports** 2D axi-symmetrical configurations and **is compatible** with `Pb_Multiphase`.

![Scheme of a VDF grid: scalars are stored at the center of the elements and the normal component of the velocities at the faces](figures/VDF-cut.png)
