@page Dev_Lata LATA Format Specification

Version 2.0 — Benoit Mathieu — 25 March 2010

---

## 1. Introduction

This document describes the LATA 2.0 post-processing format, introduced in
version 1.5.6 of Trio_U.

The LATA format is designed to record meshes and fields at successive time steps
during a computation, for post-processing purposes. Its intent is not to support
computation restarts; the stored values may therefore carry lower precision than the
internal solver representation.

The design of the LATA format is guided by the following objectives:

- Provide an extremely simple storage format, readable by a program of only a few
  lines (hence the use of an ASCII master file and minimally structured ASCII or
  binary data files).
- Support the storage of very large data volumes (hence the use of one separate
  file per field per time step).
- Offer straightforward support for dynamic meshes (use of a distinct mesh at each
  time step).
- Provide an extensible format (addition of face-located fields, face and edge
  descriptions, extra geometric element attributes, structured mesh storage, dynamic
  refinement, etc.).

Two levels of abstraction are distinguished for data written to a LATA file:

- **LataDB** (LataDataBase): a low-level container comparable to HDF, specifying
  how hierarchically organized arrays of integers or floating-point numbers are
  stored. The LataDB specification is designed to be stable across code versions.
- **Meshes and fields** (higher-level objects) are stored within the LataDB
  container according to conventions that may evolve between code versions,
  in particular through the addition of new object types.

This document is divided into two parts:

1. The specification of the LataDB format (referred to as Lata V2).
2. A description of how Trio_U uses this format to store meshes, fields, boundary
   information, and parallel computation metadata.

---

## 2. LataDB Format Definition

The reference implementation for interpreting LATA files is provided in the
`lata2dx` code (see the `LataDB` class).

Data in LataDB format is spread across multiple files:

- A **master file**, in ASCII format, containing directory information: the list of
  domains, time steps, and fields; array sizes; and the names of the data files.
- One or more **data files**, in ASCII or binary format.

### 2.1 Data File Specification

A data file contains one or more blocks, each holding an *n* × *m* matrix of
integer or floating-point values, where *n* is the number of rows and *m* is the
number of columns.

Each block is characterized by:

- The absolute byte offset of the beginning of the block within the file.
- The value encoding: ASCII, big-endian binary, or little-endian binary. With
  big-endian encoding, the most significant byte is stored at the lowest memory
  address.
- The value type: `int32`, `int64`, `ieee-real32`, or `ieee-real64`.
- The array dimensions: values *n* and *m*.
- The storage order: Fortran order (column-major, values written column by column)
  or C order (row-major, values written row by row).
- The presence and type of Fortran block markers: none, a single marker pair
  enclosing the entire matrix, or one marker pair per column.
- The binary type of Fortran block markers: 32-bit or 64-bit integer.

When present, a Fortran block marker corresponds to the value that a Fortran `READ`
statement would expect for a binary file. Specifically, it is the number of bytes of
data content. Note that when the data block is encoded as ASCII, the marker is also
written in ASCII, and its value is the number of bytes the interior of the block
*would* occupy if it had been written in binary (excluding the markers themselves).

### 2.2 Master File Specification

The LATA master file is an ASCII file that contains a hierarchical database of
pointers to data blocks located in separate data files.

#### 2.2.1 Header

The header consists of exactly three mandatory lines, optionally followed by a
fourth line specifying the default encoding format.

The first word of the first line identifies the low-level format name (`LATA_V2.1`).
The remainder of the first line, the second line, and the first word of the third
line are application-defined. For reference, Trio_U writes:

```
LATA_V2.1 Trio_U Version1.4.9
nom_du_cas
Trio_U
Format LITTLE_ENDIAN,INT32,F_INDEXING,C_ORDERING,F_MARKERS_SINGLE,REAL32
```

@note If the file does not begin with `LATA_V2.`, it is assumed to be written in
the legacy format.

The remainder of the third line is interpreted by `lata2dx` and may contain
command-line-style options. For example, the following line forces reconnection of
parallel meshes into a single block and enables verbose debug output during file
reading in VisIt:

```
Trio_U reconnect=1e-6 verbosity=4
```

The `Format` keyword defines the default encoding for all fields in the file. This
default may be overridden on a per-field basis. The recognized keywords following
`Format` are:

| Keyword | Meaning |
|---|---|
| `BIG_ENDIAN`, `LITTLE_ENDIAN`, `ASCII` | Binary byte order or ASCII encoding. Intel processors use little-endian. |
| `INT32`, `INT64` | Default size for integer fields. |
| `REAL32`, `REAL64` | Default precision for floating-point fields. Note: `REAL64` has limited support in `lata2dx`. |
| `F_INDEXING`, `C_INDEXING` | For integer fields containing indices: whether indices are 1-based (Fortran) or 0-based (C). For example, the element connectivity table contains vertex indices, numbered from 1 in Fortran and from 0 in C. |
| `F_ORDERING`, `C_ORDERING` | For multi-component arrays (vertex coordinates, vector fields, etc.): whether values are stored in C order (all components for each point in sequence) or Fortran order (all points for the first component, then all points for the second, etc.). |
| `F_MARKERS_NO`, `F_MARKERS_SINGLE`, `F_MARKERS_MULTIPLE` | Whether Fortran block markers are absent (`NO`), present as a single pair enclosing all components (`SINGLE`), or present as one pair per component (`MULTIPLE`). |

#### 2.2.2 Master File Body

The body of the master file consists of a sequence of `TEMPS`, `GEOM`, and `CHAMP`
blocks. Keywords and values may be separated by any combination of spaces, tabs, and
newlines. It is recommended to place each `TEMPS`, `GEOM`, or `CHAMP` block on a
single line to simplify parsing by scripts.

Example:

```
GEOM dom type_elem=HEXAEDRE
CHAMP SOMMETS post.lata.dom geometrie=dom size=3745182 composantes=3
CHAMP ELEMENTS post.lata.dom geometrie=dom size=3538944 composantes=8
        format=INT32 file_offset=44942192
CHAMP FACES post.lata.lata.dom geometrie=dom size=10819592 composantes=4
        format=INT32 file_offset=158188408
CHAMP ELEM_FACES post.lata.lata.dom geometrie=dom size=3538944 composantes=6
        format=INT32 file_offset=331301888
TEMPS 0
CHAMP VITESSE post.lata.lata.VITESSE.FACES.dom.pb.0.0000000000
        geometrie=dom size=10819592 composantes=1 localisation=FACES nature=scalar
CHAMP TEMPERATURE post.lata.lata.TEMPERATURE.ELEM.dom.pb.0.0000000000
        geometrie=dom size=3538944 composantes=1 localisation=ELEM nature=scalar
TEMPS 5.5923e-06
CHAMP VITESSE post.lata.lata.VITESSE.FACES.dom.pb.0.0000055923
        geometrie=dom size=10819592 composantes=1 localisation=FACES nature=scalar
CHAMP TEMPERATURE post.lata.lata.TEMPERATURE.ELEM.dom.pb.0.0000055923
        geometrie=dom size=3538944 composantes=1 localisation=ELEM nature=scalar
FIN
```

In this example, a geometry named `dom` is defined as valid for all time steps
(the `GEOM` block appears before the first `TEMPS` block), discretized with
hexahedral elements. The vertex coordinates are declared immediately after,
followed by the element connectivity (eight vertex indices per hexahedral cell,
stored as `INT32`), the face connectivity (vertex indices per face), and the
element-to-face connectivity (six face indices per element). Two time steps are
then declared, each providing a scalar normal-velocity field at faces and a
scalar temperature field at elements.

The `FIN` keyword is optional.

#### 2.2.3 "TEMPS" Entry

The `TEMPS` block declares a new time step. All geometry and field definitions that
follow belong to this time step. The physical time is specified immediately after
the keyword:

```
TEMPS 5.5923e-06
```

#### 2.2.4 "GEOM" Entry

The `GEOM` block declares a new geometry name. The current VisIt plugin supports
two cases only: static geometries, declared before the first `TEMPS` block and
valid for all time steps; and dynamic geometries, which change at every time step
and must be fully re-declared (including vertex coordinates, element connectivity,
etc.) inside each `TEMPS` block.

Syntax:

```
GEOM name type_elem=type
```

The element type string is application-defined metadata and is not interpreted by
the `LataDB` low-level container.

#### 2.2.5 "CHAMP" Entry

The `CHAMP` block is used to store all array types: vertex coordinates, element
connectivity, face connectivity, parallel joint connectivity, and physical field
values.

Syntax:

```
CHAMP field_name data_file_name
    geometrie=geometry_name
    size=n
    [ composantes=m ]
    [ localisation=SOM|ELEM|FACES ]
    [ nature=scalar|vector ]
    [ format=... ]
    [ file_offset=... ]
    [ reference=... ]
    [ noms_compo=name1,name2,... ]
```

The parameters are defined as follows:

| Parameter | Description |
|---|---|
| `field_name` | Name of the field. Certain names are interpreted specifically by the VisIt plugin (e.g., `sommets`, `elements`, `faces`, `elem_faces`, `sommets_ijk_i`). This is not part of the LataDB low-level specification. |
| `data_file_name` | Path to the file containing the array data. |
| `geometrie=name` | Name of the geometry to which this field belongs. |
| `size=n` | Number of rows in the value array. |
| `composantes=m` | Number of columns in the value array. Defaults to 1. |
| `localisation=...` | Discretization location of the field, when applicable: `SOM` (vertices), `ELEM` (elements), or `FACES`. |
| `nature=...` | Defaults to `scalar`. Set to `vector` to indicate that the field should be treated as a vector field in visualization tools; the number of components must then equal the spatial dimension. |
| `format=...` | Per-field override of the default encoding format. Uses the same keyword syntax as the header `Format` line. |
| `file_offset=N` | Byte offset of the data block from the beginning of the data file. Defaults to zero. Long integers must be used to support files larger than 2 GB. |
| `reference=...` | When the field contains indices into another field (element connectivity, face connectivity, or other), this keyword identifies the referenced field. Not currently exploited. |
| `noms_compo=...` | Comma-separated list of component names. |

---

## 3. Storage of Meshes and Fields in LataV2

This chapter specifies how Trio_U, VisIt, and `lata2dx` use the LataV2 low-level
container to store meshes and fields.

### 3.1 Structured Meshes

Structured meshes are described by three fields providing the X, Y, and Z
coordinates of the mesh vertices. The `size` attribute of each coordinate field
implicitly defines the mesh dimensions.

An optional integer field named `INVALID_CONNECTIONS`, containing 0 or 1 for each
cell, marks invalid cells. The `NO_INDEXING` format attribute must be specified for
this field to prevent the C++ plugin from subtracting one from each value (which
would occur if the default `F_INDEXING` attribute were applied).

Example for a mesh of dimensions n*i* = 192, n*j* = 192, n*k* = 96:

```
LATA_V2.1
lata2dx
Trio_U
Format LITTLE_ENDIAN,INT32,F_INDEXING,C_ORDERING,F_MARKERS_SINGLE,REAL32
GEOM dom_IJK type_elem=HEXAEDRE
CHAMP SOMMETS_IJK_I titi.lata.SOMMETS_IJK_I.dom_IJK
        geometrie=dom_IJK size=193 composantes=1
CHAMP SOMMETS_IJK_J titi.lata.SOMMETS_IJK_J.dom_IJK
        geometrie=dom_IJK size=193 composantes=1
CHAMP SOMMETS_IJK_K titi.lata.SOMMETS_IJK_K.dom_IJK
        geometrie=dom_IJK size=97 composantes=1
CHAMP INVALID_CONNECTIONS titi.lata.INVALID_CONNECTIONS.dom_IJK
        geometrie=dom_IJK size=3538944 composantes=1
        localisation=ELEM format=INT32,NO_INDEXING
TEMPS 5.5923e-06
CHAMP VITESSE titi.lata.dom_IJK_VITESSE_FACES.1
        geometrie=dom_IJK size=3613153 composantes=3
        localisation=FACES nature=vector
CHAMP TEMPERATURE titi.lata.dom_IJK_TEMPERATURE_ELEM.1
        geometrie=dom_IJK size=3538944 composantes=1
        localisation=ELEM nature=scalar
FIN
```

The supported element types for structured meshes are `HEXAEDRE` and `QUADRANGLE`.

Vertex, element, and face indices are ordered such that the *i*-index varies
fastest and the *k*-index varies slowest.

The `SOMMETS_IJK_I`, `SOMMETS_IJK_J`, and `SOMMETS_IJK_K` arrays must contain
n*i* + 1, n*j* + 1, and n*k* + 1 values respectively (Fortran block markers not
counted).

Element-located arrays must contain n*i* × n*j* × n*k* values.

Face-located arrays must have as many components as spatial dimensions (one scalar
component per face direction) and as many rows as the total number of mesh vertices.
The row index of a given face is the index of its lowest-numbered vertex; the column
index is the direction of the face normal (x, y, or z). Some row-column combinations
are therefore unused.

@note In the above example the storage order is `C_ORDERING`, so velocity
components for the x, y, and z face directions are interleaved in the data file.
If the components are stored by direction (all x-faces, then all y-faces, then all
z-faces), `F_ORDERING` must be specified instead.

### 3.2 Unstructured Meshes

The supported element types for unstructured meshes are:

- `HEXAEDRE`
- `QUADRANGLE`
- `SEGMENT`
- `TETRAEDRE`
- `TRIANGLE`
- `POLYEDRE`

For the `POLYEDRE` type, the maximum number of vertices per polyhedron equals
the number of components declared for the `ELEMENTS` field. For polyhedra with
fewer vertices than the maximum, unused vertex slots in the `ELEMENTS` array must
be filled with the value −1.

The vertex ordering within elements follows the Trio_U convention (to be detailed
in a separate document).

It is permissible to omit the element type declaration. In that case, only a point
cloud is created, and the only valid field location is `SOM`.

The following fields are recognized for mesh description. Only `SOMMETS` and
`ELEMENTS` are mandatory (`ELEMENTS` is omitted only for point clouds without a
declared element type). The `FACES` and `ELEM_FACES` fields are required if any
face-located fields exist. The joint fields are used to exploit the domain
decomposition during parallel visualization. All arrays are integer arrays except
`SOMMETS`.

| Field name | Content |
|---|---|
| `SOMMETS` | Vertex coordinates. `size` = number of vertices; `composantes` = spatial dimension. |
| `ELEMENTS` | Vertex indices of each element. `size` = number of elements; `composantes` = number of vertices per element. |
| `FACES` | Vertex indices of each face. `size` = number of faces; `composantes` = number of vertices per face. |
| `ELEM_FACES` | Element-to-face connectivity. `size` = number of elements; `composantes` = number of faces per element. |
| `JOINTS_SOMMETS` | Parallel joint data for vertices. `size` = number of MPI processes; `composantes` = 2. Column 0 contains the index of the first vertex belonging to process *i*. |
| `JOINTS_ELEMENTS` | Parallel joint data for elements. Same structure as `JOINTS_SOMMETS`. |
| `JOINTS_FACES` | Parallel joint data for faces. Same structure as `JOINTS_SOMMETS`. |

