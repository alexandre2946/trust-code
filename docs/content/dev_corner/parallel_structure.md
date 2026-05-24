@page Dev_ParallelStructure Parallel Structure

This page describes the architectural foundations of TRUST's
parallelism: the programming model, the distributed mesh data
structures, the distributed-array machinery, and the partitioning
toolkit.  It complements the user-facing @ref GG_Parallel guide —
which covers *how to run* a parallel calculation — by explaining
*how the kernel is wired* internally, for anyone touching the C++
sources.

@note This page is adapted from a 2006 design document by Benoît
Mathieu (Trio_U 1.5.1 beta).  The architectural concepts described
here are stable foundations that survived the 2015 split into TRUST
and TrioCFD and are still active today, but several class names were
renamed as part of the long migration from Trio_U → TRUST.  Notable
shifts (already incorporated into this page):
`Zone` → `Domaine`,
`Joint_Item` → `Joint_Items`,
`enum Type_Item` → `enum class JOINT_ITEM`.  Distributed arrays are
now the `TRUSTTab` family (`IntTab`, `DoubleTab`, and their `_t`
templated forms `IntTab_t`/`DoubleTab_t` for 32/64-bit integer
support).  When a specific identifier matters in your work, `grep`
in `$TRUST_ROOT/src/Kernel/` to confirm the current spelling.

## Programming model: SPMD with fine-grained MPMD blocks

TRUST runs in **SPMD** mode (Single Program, Multiple Data): every
processor executes the same code path, and at any instant they sit at
the same point in the code.  All processors therefore solve the same
system of equations simultaneously (same equations, same
discretization, same time scheme).  Whenever processors need to
exchange data, they do so at **synchronisation points** — collective
operations baked into the code.

SPMD contrasts with MPMD (Multiple Program, Multiple Data), where each
processor runs a different code and inter-processor communications go
through *rendez-vous points* — less frequent but harder to organise.
For example, a hypothetical MPMD TRUST could have one processor solve
Navier–Stokes while another solves a thermal problem.

In practice, TRUST is **not** strictly SPMD: the data each processor
operates on is not always the same size (one processor's sub-domain
has a different number of elements than another's), and some
non-linear models (turbulence, convection operators, ...) have
conditional branches whose outcome depends on local values and can
therefore differ across processors.

TRUST resolves this with **fine-grained MPMD blocks** — short
sections (loops over geometric items whose count varies per
processor, loops over matrix entries, ...) that *cannot* contain
inter-processor communications.  These blocks are themselves
organised in an SPMD pattern: every processor enters the same MPMD
block at the same point in the code.

### Developer convention

The TRUST kernel's high-level algorithms (operators, matrix solvers,
framework methods) are coded SPMD by default.  When you write such
code:

- if you call a service inside a conditional (`if`, `while`, `for`),
  the conditional must produce **the same outcome on every
  processor** — same branch taken, same number of iterations.  A
  divergence here will deadlock the next synchronisation point.
- a service is implicitly SPMD as soon as it contains synchronisation
  points.  This property is not currently documented in the code —
  it's on the developer to verify by reading the implementation.
- a block of code whose conditional branches *can* legitimately
  differ per processor must be marked as MPMD (no synchronisation
  points inside it).

### Process groups (mostly unused today)

TRUST theoretically supports a higher-level MPMD layer via the
`Process` class and the notion of **groups** (sets of processors
executing SPMD code together).  A group can be the audience of an
`mpsum` call, for instance, so the sum spans only that group's
processors.  The mechanism is plumbed but has only ever been used in
prototype form (an AMR/dynamic-load-balancing experiment).  Routinely
using it would require generalising every distributed object to
carry a reference to its owning group.

## Distributed mesh and geometric items

The computational domain is partitioned into N sub-domains: every
element belongs to exactly one sub-domain.  Each processor owns one
sub-domain and knows about the geometric items attached to it
(vertices, faces, edges, boundary faces).  Those are the
sub-domain's **real items**.

Some real items are shared with a neighbouring sub-domain — a vertex
on the boundary between two sub-domains is real on both.  Such
shared items are called **common items** ("items communs").  The
boundary between two sub-domains — the set of common items plus the
metadata describing the link — is called a **joint**.

For numerical schemes computed on items adjacent to joints, a
processor needs to access the values on items just *across* the
joint, on the neighbouring sub-domain.  Those extra items are
**virtual items** ("items virtuels").  The number of layers of
virtual items kept is the **joint thickness** ("épaisseur de joint"):

- thickness 1 covers 1st- and 2nd-order space discretizations;
- thickness 2 covers 3rd- and 4th-order schemes.

If the thickness is N, a sub-domain knows every element whose
neighbours-of-neighbours-of-... (N-1 hops) include one of its real
elements, plus the attached vertices/faces/edges.  Two elements are
"neighbours" if they share at least one common vertex.

Two conventions worth knowing when reading the code:

- common items are *also* real items;
- boundary faces ("faces de bord") are *also* faces.

## Distributed arrays

A **distributed array** ("tableau distribué") holds values indexed
by geometric items.  For values associated with virtual items, the
array stores a copy obtained by copying the corresponding real value
from the neighbour processor's sub-domain.

The conventional layout is:

- **real values first**, in the order in which the geometric items
  themselves are defined (entry `n` of the array is the temperature
  of element `n` in the local sub-domain, say);
- **virtual values last**, ordered first by the index of the
  sub-domain that owns the corresponding real item, then within
  each sub-domain by the order defined by the **distant items**
  table (the list of *its* real items that a neighbour is mailing
  out).

### Virtual space exchange

The standard pattern: every processor computes the values on its
real items, then performs an **"échange espace virtuel"** (virtual
space exchange) where it sends those real values to the neighbours
that need them as virtual values, and receives the values its own
virtual items expect.  After this round-trip, the array's virtual
section is up to date.

A few invariants:

- a real item can be sent to several processors but is always
  *received* from exactly one;
- a single item can be both *distant* (sent to others) and *common*
  (shared with others as a real);
- for **common items** specifically, every co-owner could in
  principle recompute the same value.  In exact arithmetic, with
  identical inputs, this would converge.  In practice rounding
  errors break that — so TRUST performs an additional **"échange
  des items communs"** in which the **lowest-rank sub-domain that
  owns the item imposes its value** on every other owner.  This
  keeps the parallel solution bit-identical across processors.

### Constructing distributed arrays

The virtual-space layout of a distributed array is described by an
`MD_Vector` (metadata vector) — built once per geometric type for a
domain and reused by every array indexed by items of that type.  The
`Scatter` class exposes the construction sequence:

- `Scatter::calculer_espace_distant_*` (one variant per geometric
  type — `_elements`, `_sommets`, `_faces`, `_aretes`) computes, for
  this sub-domain, the list of real items it must export as virtuals
  to each neighbour;
- `Scatter::construire_md_vector(domaine, nb_items_reels, JOINT_ITEM, md)`
  packages the per-neighbour common + distant lists plus the virtual
  count into an `MD_Vector`;
- the array is then created with that `MD_Vector` attached — every
  subsequent `echange_espace_virtuel()` call uses it to drive the
  MPI exchange.

  ![Example of array parallel structure](figures/parallel_example.png)

  
### Arrays of geometric indices: index translation

Many distributed arrays hold not numerical values but *indices into
other geometric arrays* — e.g. for each element, the local indices
of its vertices in the local-sub-domain vertex table.  Naïvely
running the virtual-space exchange on such an array would copy raw
indices from another processor, which would point at the *wrong*
vertices on the receiving side.

`Scatter::construire_espace_virtuel_traduction(md_indice, md_valeur,
tableau)` handles this: it exchanges the raw indices via the
`md_indice` layout, then translates them into the local numbering
defined by `md_valeur`.  The element-vertex table, the face-vertex
table, the element-faces and face-neighbours tables all rely on
this mechanism.

## The "Joint" and "Joint_Items" data structures

Each sub-domain stores its joints in a list of `Joint` objects, one
per neighbouring sub-domain.  Each `Joint` then carries one
`Joint_Items` instance per geometric type, indexed by the
`enum class JOINT_ITEM { SOMMET, ELEMENT, FACE, ARETE, FACE_FRONT }`
declared in `Joint.h`.  A `Joint_Items` exposes:

| Accessor                 | What it returns |
|--------------------------|---|
| `items_distants()`       | Local indices of the real items this sub-domain exports as virtuals to the neighbour (`ArrOfInt_t`). |
| `items_communs()`        | Local indices of the items shared (real on both sides) with the neighbour (`ArrOfInt_t`). |
| `renum_items_communs()`  | Two-column index map (`IntTab_t`): column 0 — local index; column 1 — corresponding index on the neighbour processor. |

For initialisation, the matching `set_items_distants()`,
`set_items_communs()` and `set_renum_items_communs()` setters return
mutable references.  All accessors come in templated `_32_64`
flavours — the `_t` typedefs (`ArrOfInt_t`, `IntTab_t`) resolve to
the right integer width at compile time.

## Distributed geometry specifications

For each canonical mesh array, here's the convention:

- **Vertex coordinates** — `Domaine::les_sommets()` returns the
  `DoubleTab_t&` of 2D/3D vertex coordinates; `coord_sommets()` is
  the const accessor.  Standard distributed array with common items
  + virtual space, populated by the standard virtual-space exchange.

- **Element table** — `Domaine::les_elems()` returns an `IntTab_t&`
  giving, for each cell, the local indices of its vertices.
  Standard distributed array with virtual space but *no common
  items*; virtual contents obtained via the index-translation
  exchange.

- **Face-vertex table** — `Domaine_VF::face_sommets()` returns an
  `IntTab&` giving, per face, its local vertex indices.  Vertex
  order within each face is determined by
  `Elem_geom_base::get_tab_faces_sommets_locaux()`.  Standard
  distributed array with virtual space + common items; virtual
  contents obtained via index translation.

- **Element-face table** — `Domaine_VF::elem_faces()` returns an
  `IntTab&` giving, per cell, its face indices.  Standard distributed
  array indexed by elements with virtual space.

- **Face neighbours table** — `Domaine_VF::face_voisins()` returns
  an `IntTab&` giving, per face, the indices of its neighbouring
  elements (`-1` when the face has only one neighbour).  Standard
  distributed array indexed by faces with virtual space + common
  items.  Discretisation-specific conventions apply once exchanged:
  - **VDF**: `face_voisins(i,0)` is the element with the smaller
    coordinate in the face's normal direction; `face_voisins(i,1)`
    is the larger.  For boundary faces and virtual faces, one of
    the two neighbours is `-1`; the order of `voisin0` / `voisin1`
    is identical on every processor that knows the face.
  - **VEF**: no ordering constraint between `face_voisins(i,0)` and
    `face_voisins(i,1)` for faces with two neighbours.  For
    one-neighbour faces, the known neighbour is in
    `face_voisins(i,0)` and the other slot is `-1`.  ⚠️ the order
    of neighbours is *not* guaranteed identical across processors
    sharing the same face.

## The partitioning toolkit

The mesh cutter is organised around three classes:

- `Decouper` — an `Interprete` that reads the partitioner's
  parameters from the data file and drives the operations
  (partition the domain, write the `.Zones` files, write the
  partitioning file, ...).
- `Partitionneur_base` — abstract base of the partitioner
  hierarchy.  Concrete subclasses implement the individual
  algorithms (Metis, Tranche, ...).
- `DomaineCutter` — given a complete domain plus a partition map,
  builds the per-sub-domain data: real vertices, vertex/element
  tables with local indices, sub-domain borders, joint vertices,
  joint faces.  (It also computes the distant items for now, though
  that step is optional and can be deferred to `Scatter` — useful
  in particular when reading a parallel MED file that already
  contains the partition but not the distant items.)

### Constraints on every partitioner

`Partitionneur_base::construire_partition()` must produce a
partition that respects periodic boundaries.  Two invariants:

- **If a real periodic face is on a processor, the opposite
  periodic face must also be on that processor** — i.e. the
  element across the periodic boundary lands in the same
  sub-domain.
- **The `renum_som_perio` table must point at the same real vertex
  on every processor.**  If a vertex is multi-periodic, every
  processor that owns it must also own every other multi-periodic
  vertex associated with it.  A partition that splits a
  multi-periodic vertex group across processors is rejected.

  Similarly, **a processor must not own an isolated periodic
  vertex**: if a sub-domain contains a vertex adjacent to a
  periodic face, it must own at least one of the adjacent periodic
  faces, so the periodic vertex can be discovered by walking the
  faces.

To enforce these, `Partitionneur_base` exposes
`corriger_bords_avec_graphe()` and `corriger_bords_avec_liste()`,
which post-process a candidate partition until the invariants hold:

- if a sub-domain has a periodic face, it also has the
  corresponding opposite face;
- no isolated periodic vertex on any sub-domain;
- multi-periodic vertex groups stay together.

### What "DomaineCutter" does per sub-domain

For each sub-domain it builds:

- the set of real vertices (every vertex of every element assigned
  to the sub-domain);
- the vertex table and element table with local indices;
- the sub-domain's borders;
- the **joint vertices** (vertices shared with each neighbouring
  sub-domain), ordered so they share the same order on the
  neighbouring side;
- the joint faces.

The distant items are also computed here in the current code path,
but that operation is optional — it can be deferred to `Scatter` at
load time (notably for MED-format partitions that don't ship the
distant-element lists).
