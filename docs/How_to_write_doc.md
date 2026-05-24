# How to write documentation for TRUST

This guide is for contributors who want to add or edit pages in the TRUST
documentation. It assumes you know basic Markdown but have not used Doxygen
before.

For build instructions, directory layout, and a map of every file in the
system, see `README.md` in this same directory.

---

## Table of contents

1. [Navigation tree](#1-navigation-tree)
2. [Standard Markdown](#2-standard-markdown)
3. [Math formulas](#3-math-formulas)
4. [Figures](#4-figures)
5. [Cross-page links](#5-cross-page-links)
6. [Callout boxes](#6-callout-boxes)
7. [Bibliography and citations](#7-bibliography-and-citations)
8. [Checklist for a new page](#8-checklist-for-a-new-page)

---

## 1. Navigation tree

### Page declaration

Every `.md` file that should appear as a page in the documentation **must**
begin with a `@page` declaration:

```
@page MyPageID  Title shown in the navtree
```

- **`MyPageID`** — a unique identifier used to reference this page elsewhere.
  Use `CamelCase` and pick something descriptive (e.g. `Dev_MyTopic`,
  `Tut_MyTutorial`).
- **Title** — free text, shown in the sidebar and browser tab.

### Registering a new page in the tree

A page not referenced by any other page becomes an orphan: it is built but
not reachable from the sidebar. To attach it, open the appropriate parent
index file and add a `@subpage` line:

```markdown
@page SomeParentIndex  Parent section title

- @subpage MyPageID
- @subpage AnotherPageID
```

`@subpage` both creates a hyperlink **and** declares the page as a child in
the navtree — it does both jobs at once.

### Where to add your page

For the full map of what goes where (narrative pages, keyword descriptions,
trustutils API, sidebar config, etc.) see the **Adding or modifying
documentation** table at the bottom of `README.md`.

### Example: adding a new developer page

1. Create `doc_md/dev_corner/my_topic.md`:

```markdown
@page Dev_MyTopic My topic title

Content goes here ...
```

2. Open `doc_md/dev_corner/index.md` and add:

```markdown
- @subpage Dev_MyTopic
```

That is all that is required.

---

## 2. Standard Markdown

All standard Markdown features work as expected.

### Headings

```markdown
## Section
### Subsection
#### Sub-subsection
```

> **Do not use `# H1`** inside a page body — Doxygen reserves the top-level
> heading for the page title set in `@page`. Start at `##`.

> **Avoid `## headings` on the mainpage alongside `@subpage` links.**
> Doxygen registers every Markdown `##` heading as a navtree entry.  On the
> mainpage (`main_page.md`) this causes the section to appear in the sidebar
> — typically nested under the first `@subpage` child, which is wrong.
> Use a raw HTML tag instead to get the visual heading without the navtree
> entry:
>
> ```html
> <h2>My Section</h2>
> ```

### Emphasis, inline code

```markdown
*italic*   **bold**   `inline code`
```

### Lists

```markdown
- item one
- item two
  - nested item

1. first
2. second
```

### Tables

```markdown
| Column A | Column B |
|---|---|
| value    | value    |
```

### Code blocks

Fenced blocks with an optional language tag for syntax highlighting:

````markdown
```cpp
double x = 1.0;
```
````

Supported tags include `cpp`, `python`, `bash`, `markdown`.

---

## 3. Math formulas

TRUST documentation uses **MathJax 3** (LaTeX syntax). Doxygen uses its own
delimiters rather than `$...$` to avoid clashing with dollar signs in code
examples.

### Inline math

Wrap with `\f$` ... `\f$`:

```
The velocity field \f$\boldsymbol{u}\f$ satisfies ...
```

### Display (block) math

Wrap with `\f[` ... `\f]` on their own lines:

```
\f[
\partial_t \boldsymbol{u} + \nabla \cdot \boldsymbol{F} = \boldsymbol{S}
\f]
```

### Multi-line aligned equations

Use the standard LaTeX `aligned` environment inside `\f[` ... `\f]`:

```
\f[
\left\{
\begin{aligned}
a_h(\boldsymbol{u}_h, \boldsymbol{v}_h) + b_h(\boldsymbol{v}_h, p_h) &= L_h(\boldsymbol{v}_h) \\
c_h(\boldsymbol{u}_h, q_h) &= 0
\end{aligned}
\right.
\f]
```

### Common LaTeX snippets

| Construct | LaTeX |
|---|---|
| Bold vector | `\boldsymbol{u}` |
| Partial derivative | `\partial_t` |
| Divergence / gradient | `\nabla \cdot`, `\nabla` |
| Integral | `\int_{\omega}` |
| Fraction | `\frac{a}{b}` |
| Subscript / superscript | `u_f^{n+1}` |
| Norm | `\|\boldsymbol{u}\|` |

---

## 4. Figures

### Storing images

Place all image files in `docs/content/figures/`. The Doxyfile sets
`IMAGE_PATH = figures`, so Doxygen finds images by bare filename regardless
of which subdirectory the `.md` source file is in.

Preferred format: **PNG**. JPEG and SVG are also supported.

### Inserting a figure

```markdown
![Brief description of the image](figures/my_image.png)
```

The alt text is rendered as a caption below the image. Example from the VEF
page:

```markdown
![Control volume for velocity](figures/control_volume_velocity.png)
```

### Controlling the size

Doxygen does not support Markdown image sizing. If you need to control the
rendered width, use an HTML `<img>` tag:

```html
<img src="figures/my_image.png" alt="Description" width="400"/>
```

---

## 5. Cross-page links

### Hyperlink to another page

```
See the @ref DevCorner "Developer Corner" for more details.
```

The quoted text is the link label. If omitted, Doxygen uses the page title.

### Link to a section on another page

Doxygen auto-generates anchors from Markdown headings (lowercase, spaces
replaced by underscores):

```
See @ref Disc_VEF#finite_volume_approach "the VEF finite volume approach".
```

### Referencing a C++ class

Doxygen auto-links known C++ identifiers when written in backticks. For an
explicit link with custom label:

```
@ref Equation_base "base equation class"
```

---

## 6. Callout boxes

Doxygen provides special-paragraph commands that render as styled boxes.
Place them on their own paragraph (blank line before and after).

| Command | Rendered style | Typical use |
|---|---|---|
| `@note` | Blue info box | Important clarifications |
| `@warning` | Orange warning box | Pitfalls, breaking changes |
| `@todo` | Purple todo box | Known gaps, future work |
| `@attention` | Yellow attention box | Things the reader must not miss |

Example:

```
@note In the Cavite1 block, \f$L_x = 0.3\f$ and the number of nodes is 4
along X.

@warning TRUST is sensitive to whitespace. Use a space before and after
each keyword.
```

---

## 7. Bibliography and citations

### Citing a reference

Use `@cite key` anywhere in the text, where `key` is the BibTeX entry
identifier in `docs/content/references.bib`:

```
Initially introduced in @cite LM89, the VEF method ...
```

Multiple consecutive citations are written one after the other:

```
This was studied in @cite Heib2003 @cite Fortin2006.
```

Doxygen renders each `@cite` as a numbered superscript linked to the
Bibliography page, and builds the full reference list automatically.

### Adding a new reference

Open `docs/content/references.bib` and append a standard BibTeX entry:

```bibtex
@article{AuthorYYYY,
  author  = {Last, First and Other, Name},
  title   = {Full title of the article},
  journal = {Journal Name},
  year    = {2024},
  volume  = {42},
  pages   = {1--20},
  doi     = {10.xxxx/xxxxx},
}
```

Common entry types: `@article`, `@book`, `@inproceedings`, `@techreport`,
`@misc`.

The key convention used in this project is `FirstAuthorLastNameYY`
(two-digit year), e.g. `LM89`, `Heib2003`.

### How the Bibliography navtree entry works

Doxygen auto-generates the bibliography as `citelist.html` but **never puts
it in the navtree** — even `@subpage citelist` is silently ignored.  The
workaround is a two-part mechanism:

1. `doc_md/references.md` declares `@page References Bibliography`.  This is
   a normal Doxygen page that Doxygen *does* wire into the navtree.  It must
   not be deleted or renamed.
2. `main_page.md` lists `@subpage References`, placing it at the top level.
3. `flatten_navtree.py` post-processes the generated JS to replace every
   `"References.html"` URL in the navtree data with `"citelist.html"`, so
   clicking "Bibliography" in the sidebar lands directly on the full reference
   list rather than the placeholder page.

**Do not add content to `references.md`.**  It exists only as a navtree
placeholder; users never see it directly.

---

## 8. Checklist for a new page

Before pushing, verify:

- [ ] The file starts with `@page UniqueID Title`.
- [ ] The page ID is referenced with `@subpage` in the appropriate parent
      index file.
- [ ] The file is **not** listed in `EXCLUDE_PATTERNS` in `Doxyfile` (check
      if `References.html` not appearing in the build is the symptom).
- [ ] Figures are stored in `docs/content/figures/` and referenced as
      `figures/name.png`.
- [ ] Math uses `\f$...\f$` (inline) or `\f[...\f]` (block) — **not** `$...$`.
- [ ] New bibliography entries are added to `docs/content/references.bib`.
- [ ] The build completes cleanly: run `make all` from `docs/`
      and check the terminal output for warnings about undefined references
      or missing images.
- [ ] Check both light and dark mode visually — if a new element appears with
      a hardcoded white background in dark mode, add a CSS variable override
      in `docs/theme/custom.css` (see `README.md` → Theming and dark mode).
