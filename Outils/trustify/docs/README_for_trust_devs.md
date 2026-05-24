@page Trustify_ForTrustDevs trustify — guide for TRUST/baltik developers

## Documenting your keywords with `// XD` tags

TRUST reads its `.data` input files using keywords (one per readable C++
class). **trustify** turns the `// XD` comments you write next to those
classes into a machine-readable schema that powers dataset validation,
auto-formatting, the reference manual, and editor auto-completion (LSP).

Adding the tags is cheap and pays off immediately: as soon as your
keyword is decorated, `trustify check` validates every dataset that uses
it, the LSP suggests it (and its attributes) while users type, and it
appears in the generated documentation. **If you don't decorate a new
keyword, trustify cannot see it** — datasets that use it will fail to
validate.

This page is the reference for *writing* those tags. It assumes no prior
trustify knowledge.

---

### TL;DR — a complete example

You have a class `the_kw`, deriving from an existing keyword `base_kw`,
that reads one floating-point parameter `tutu` between braces:

```cpp
// XD the_kw base_kw the_kw BRACE A one-line description of the keyword.
Implemente_instanciable(The_kw, "the_kw", Base_kw);
// ...
void The_kw::set_param(Param& param)
{
  param.ajouter("tutu", &tutu_);  // XD_ADD_P floattant A description for the tutu attribute.
}
```

That is all. Re-generate the schema (`trustify generate_schema`, or just
run `trustify check some_dataset.data`, which regenerates on demand) and
`the_kw { tutu 3.14 }` becomes valid.

The rest of this page explains every field, every attribute type, and
the handful of extra tags for lists, dictionaries, long lines, and
nested parameters.

---

### Where the tags live

The scanner only reads **`.cpp` and `.xd` files**. Tags placed in a
header (`.h`) file are **silently ignored** — put them in the `.cpp`.

A `// XD` tag is an ordinary C++ comment; it has no effect on
compilation. Two forms exist:

- inline on a `param.ajouter(...)` line (the `XD_ADD_P` form — most
  common), and
- standalone, on its own comment line (the `XD attr` form).

You can also declare keywords that have no C++ class at all in a
standalone **`.xd` file** (same grammar; see
[Keywords without a C++ class](#keywords-without-a-c-class)).

A note on case: every structural field (names, parent, flags, types) is
lowercased by the scanner, so `BRACE`/`brace` and `The_Kw`/`the_kw` are
equivalent. Descriptions keep their original case.

---

### 1. Declaring the keyword — the header line

```
// XD <name> <name_base> <syno> <brace_flag> <description...>
```

| Field | Meaning |
|---|---|
| `<name>` | The keyword as typed in datasets, lowercase. Allowed characters: `[A-Za-z0-9_]`. |
| `<name_base>` | The **parent** keyword. Must be a keyword that is already declared (anywhere in the scanned sources) or the universal root `objet_u`. A keyword cannot be its own parent. The new keyword inherits all of the parent's attributes. |
| `<syno>` | Pipe-separated synonyms, e.g. `lire\|read`. **When there are no synonyms, repeat the keyword name** (`the_kw`). Some internal helper blocks use the literal `nul` instead — both are fine, they just mean "no real alias". |
| `<brace_flag>` | `BRACE`, `NO_BRACE`, or `INHERITS_BRACE` — see [Brace flags](#brace-flags). |
| `<description...>` | Free text. Used in error messages and the generated manual. Write a real sentence; the common placeholder `not_set` renders as an empty description. |

Synonyms registered through the `Implemente_instanciable(Class,
"syno1|syno2", Parent)` macro and `Class::add_synonym("other")` calls are
picked up automatically and merged into the `<syno>` field — you do not
need to repeat them in the `// XD` line.

---

### 2. Declaring attributes

An attribute is one parameter the keyword reads. There are two ways to
declare one.

#### 2a. Inline — `XD_ADD_P` (recommended)

Place the tag on the **same line** as the `param.ajouter(...)` call. The
attribute *name* is read from the first argument of `ajouter`, so you
only write the type and description:

```cpp
param.ajouter("tutu", &tutu_);  // XD_ADD_P floattant Description of tutu.
```

```
// XD_ADD_P <type> <description...>
```

The same applies to the whole `ajouter*` family (`ajouter_flag`,
`ajouter_non_std`, `ajouter_condition`, ...). The call and its tag
**must be on one physical line** (a multi-line `ajouter(...)` carrying an
`XD_ADD_P` is rejected) — though the *description* may overflow onto
`XD_CONT` lines (see [Long lines](#long-lines)).

**Optionality is inferred automatically**: if the `ajouter` call passes
`Param::REQUIRED`, the attribute is required; otherwise it is optional.

```cpp
param.ajouter("haspi", &haspi_, Param::REQUIRED);  // XD_ADD_P floattant Required suction height.
```

#### 2b. Standalone — `XD attr`

When there is no `param.ajouter` line to hang the tag on (hand-written
parsing, or a `.xd` file), declare the attribute on its own line,
**directly under** the block it belongs to. Here you spell out the name
and the optional/required flag yourself:

```
// XD attr <name> <type> <syno> <REQ|OPT> <description...>
```

```cpp
// XD attr iter_min entier iter_min OPT Minimum number of iterations (default 1).
```

- `<syno>` — same convention as the header line (pipe-separated; repeat
  the name for none).
- `<REQ|OPT>` — `REQ` for mandatory, `OPT` for optional. (You write this
  explicitly here; the `XD_ADD_P` form derives it from
  `Param::REQUIRED`.)

Ordering matters: an `XD attr` line attaches to whichever block header
was most recently opened in the file.

---

### 3. Attribute types

This is the part most worth getting right. The `<type>` slot accepts far
more than plain scalars.

#### Scalars

| Type | C++ shorthand | Reads in the dataset | Pydantic |
|---|---|---|---|
| `entier` | `int` | a single integer | `int` |
| `floattant` | `double` | a single float | `float` |
| `chaine` | | a single word **or** a whole brace-block (see below) | `str` |
| `rien` | `flag` | nothing — the attribute is a flag (present/absent) | `bool` |

`chaine` is deliberately permissive: `toto` is a `chaine`, and so is the
entire block `{ how do you { do } }` — but `titi toto` (two words) is
not. Use it for a single token or for an opaque brace-delimited payload.

`rien` (alias `flag`) is for switches where only the presence of the
keyword matters, e.g. `param.ajouter_flag("gas_turb", &gas_turb_); //
XD_ADD_P flag`.

#### Composition — another keyword as the type

Any declared keyword name can be used as an attribute type. This nests
one keyword inside another, which is the normal way to build up the
schema:

```cpp
param.ajouter("diffusion", &diff_);  // XD_ADD_P field_base The diffusion field.
```

Here `field_base` is itself an XD-declared keyword; the user supplies any
of its subclasses, fully parsed and validated.

#### References — `ref_<keyword>`

`ref_<keyword>` means "the *name* of an object of type `<keyword>`
declared elsewhere in the dataset", not an inline definition. It parses
as a bare string but documents the intended target:

```cpp
// XD attr domain_name ref_domaine domain_name REQ Name of the domain.
```

#### Constrained scalars

| Type | Meaning |
|---|---|
| `chaine(into=["a","b","c"])` | string restricted to the listed values (an enum) |
| `entier(into=[1,2,4])` | integer restricted to the listed values |
| `entier(min=0)` | integer `>= 0` |
| `entier(max=10)` | integer `<= 10` |
| `entier(min=0,max=10)` | integer in the closed range `[0, 10]` (`min` always precedes `max`) |

Each of these accepts an optional, **documentation-only** default after
the constraint: `chaine(into=["a","b"],default="b")` or
`entier(min=0,max=10,default=3)`. The default is shown in the manual but
does **not** change parsing behaviour (an unset optional attribute is
still `None`). The declared default must itself be a legal value, or
schema generation errors out.

Nothing may trail the closing `)` — `chaine(into=["a"])toto` is rejected.

#### Lists

Count-prefixed lists — the dataset gives the length first, then the
values (`3 1.0 2.0 3.0`):

| Type | Element |
|---|---|
| `list` | floats |
| `listentier` | integers |
| `listchaine` | strings (single words) |

Dimension-sized lists — the length is the problem `dimension` (2 or 3),
so **no count** is written in the dataset (`1.0 2.0 3.0`):

| Type | Element |
|---|---|
| `listf` | floats |
| `listentierf` | integers |
| `listchainef` | strings |

For a list **of full objects** (e.g. a list of source terms or
equations), use a `listobj` block — see
[Lists of objects](#lists-of-objects).

#### `dico` — a dictionary/enum built from C++

A `dico` is an enum whose allowed values are scraped from the C++
`dictionnaire(...)` calls. Declare the attribute `dico`, then tag each
`dictionnaire` line with `XD_ADD_DICO`:

```cpp
param.ajouter("reorder", &reorder_);   // XD_ADD_P dico How to reorder the mesh.
param.dictionnaire("none", 0);         // XD_ADD_DICO No reordering (the default).
param.dictionnaire("morton", 1);       // XD_ADD_DICO Morton-curve reordering.
param.dictionnaire("hilbert", 2);      // XD_ADD_DICO Hilbert-curve reordering.
```

This produces `chaine(into=["none","morton","hilbert"])`. The
`XD_ADD_DICO` lines must follow the `XD_ADD_P dico` line; the first
argument of each `dictionnaire(...)` becomes an allowed value.

#### `suppress_param` — drop an inherited attribute

When your keyword inherits from a parent but should *not* expose one of
the parent's attributes, suppress it:

```cpp
// XD attr gravite suppress_param gravite OPT del
```

The named attribute must actually exist in the parent, or schema
generation errors out.

---

### Brace flags

The brace flag on the header line says how the keyword is delimited in
the dataset:

- **`BRACE`** — the keyword is followed by `{ ... }` and attribute names
  are written explicitly. Example: `lire sch { tinit 0. tmax 0. }`.
- **`NO_BRACE`** — no braces; attribute *values* are read positionally,
  in declaration order, and names are not written. Example:
  `champ_uniforme 3 0. 0. 0.`.
- **`INHERITS_BRACE`** — use whatever the parent class uses. Convenient
  for a subclass that does not change the syntax.

Legacy numeric flags (`1`/`0`/`-1`, and the historical `-3`/`-2`) still
parse but trigger a modernization warning — run `trustify modernize` to
rewrite them to the named forms.

---

### Long lines

XD lines longer than 120 characters trigger a warning. Continue any XD
line with `// XD_CONT`, which must come **immediately after** the line it
continues. It works for header lines, `XD attr` lines, and `XD_ADD_P`
descriptions:

```cpp
// XD the_kw base_kw the_kw BRACE A really long and verbose description that
// XD_CONT continues here because one line was not enough, and
// XD_CONT even spills onto a third.

param.ajouter("tutu", &tutu_);  // XD_ADD_P floattant A long description for tutu that
                                // XD_CONT also needs several lines.
```

`trustify modernize` can also re-flow over-long XD lines into `XD_CONT`
continuations automatically.

---

### Nested parameters — `2XD` / `3XD`

When a keyword opens a *secondary* parameter block in C++ (via
`ajouter_param(...)` / a `Read_Sub_Param`-style helper), declare the
inner block's attributes with the `2XD` family — same grammar, one level
deeper. `3XD` goes one level deeper still. **Three levels maximum.**

```cpp
Param param("lire_ceg");
param.ajouter("t_deb", &t_deb_, Param::REQUIRED);  // XD_ADD_P floattant Initial calculation time.
param.ajouter("dt_post", &dt_post_);               // XD_ADD_P floattant Printing period, in seconds.

Param& sub = param.ajouter_param("CEA_JAEA");      // XD_ADD_P ceg_cea_jaea CEA_JAEA's criterion.
// 2XD ceg_cea_jaea objet_lecture nul NO_BRACE Parameters for the CEA_JAEA criterion.
sub.ajouter("normalise", &normalise_);             // 2XD_ADD_P entier Renormalize (1) or not (0).
sub.ajouter("nb_mailles_mini", &nb_mailles_mini_); // 2XD_ADD_P entier Minimum cells to detect a vortex.
```

The nested block (`ceg_cea_jaea` above) is declared with a `2XD` header
and its attributes with `2XD_ADD_P`. Use `2XD_CONT` / `3XD_CONT` to
continue nested lines. The scanner keeps a separate "current block" per
level, so the nested attributes attach to the right parent regardless of
interleaving.

---

### Lists of objects

A `listobj` block declares a keyword that *is* a list of other objects
(reactions, equations, sub-mediums, ...). The header line takes two
extra fields — the item type and a comma flag — and carries **no
attributes**:

```
// XD <name> listobj <syno> <brace_flag> <itemtype> <comma_flag> <description...>
```

```cpp
// XD reactions listobj nul BRACE reaction COMMA List of reactions.
```

- `<itemtype>` — the keyword type of each element (here `reaction`).
- `<comma_flag>` — `COMMA`, `NO_COMMA`, or `INHERITS_COMMA`, depending on
  whether the dataset separates elements with commas.

For lists of *builtin* scalars (floats/ints/strings) you do not need a
`listobj` block — use the `list` / `listentier` / `listchaine` attribute
types above.

---

### Keywords without a C++ class

Some keywords are pure grammar fragments with no backing class (comment
markers, fixed-size word tuples, opaque readers). Declare these in a
standalone **`.xd`** file under `src/`. The grammar is identical, except
everything uses the `// XD` / `// XD attr` standalone forms (there is no
`param.ajouter` to attach to):

```
// XD deuxmots objet_lecture nul NO_BRACE Two words.
// XD attr mot_1 chaine mot_1 REQ First word.
// XD attr mot_2 chaine mot_2 REQ Second word.
```

(See `src/Kernel/Utilitaires/generic.xd` for the canonical examples.)

---

### Checking your work

After adding or changing tags:

1. **Regenerate and surface errors.** `trustify generate_schema` (or any
   `trustify check`/`format` command, which regenerates on demand) parses
   every tag. Malformed tags are reported with their exact source
   location — and *all* errors in a pass are reported together, not just
   the first. It also warns about long lines and legacy numeric flags.
2. **Validate real datasets.** `trustify check my_dataset.data` confirms
   your keyword parses as intended; `trustify batch-check
   $TRUST_ROOT/tests` runs it across the whole test corpus.
3. **Modernize.** `trustify modernize` (dry-run by default; `--apply` to
   write) upgrades legacy flags and re-flows long lines in the current
   project's sources.

Common pitfalls the scanner will catch for you:

- A type that is neither a builtin nor a declared keyword → *"unresolved
  type ... did you forget to define it?"*
- A parent (`name_base`) that was never declared.
- `XD_ADD_DICO` without a preceding `XD_ADD_P dico`.
- `XD_CONT` not immediately following the line it continues.
- An `XD_ADD_P` on a multi-line `ajouter(...)` call (keep it on one
  line).

And one it cannot catch: **tags in `.h` files are ignored** — keep them
in `.cpp` / `.xd`.

---

### Cheat sheet

```cpp
// Header line (one per keyword):
// XD <name> <parent> <syno|name> <BRACE|NO_BRACE|INHERITS_BRACE> <desc>

// Inline attribute (name + optionality come from the C++ call):
param.ajouter("p", &p_);                 // XD_ADD_P <type> <desc>
param.ajouter("p", &p_, Param::REQUIRED);// XD_ADD_P <type> <desc>   (required)

// Standalone attribute (spell out name + REQ/OPT):
// XD attr <name> <type> <syno|name> <REQ|OPT> <desc>

// Continue a long line (immediately after):
// XD_CONT <more text>

// Dictionary enum:
param.ajouter("k", &k_);   // XD_ADD_P dico <desc>
param.dictionnaire("a", 0);// XD_ADD_DICO <desc>

// Nested blocks: prefix tags with 2 or 3 (2XD, 2XD_ADD_P, 2XD_CONT, ...)
```

Types: `entier` (`int`) · `floattant` (`double`) · `chaine` · `rien`
(`flag`) · `<keyword>` · `ref_<keyword>` · `chaine(into=[...])` ·
`entier(into=[...]|min=|max=)` · `list` · `listentier` · `listchaine` ·
`listf` · `listentierf` · `listchainef` · `dico` · `suppress_param`.
