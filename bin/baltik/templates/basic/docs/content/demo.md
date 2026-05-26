@page @@PROJECT_NAME@@_demo Getting started

This is a placeholder page for **@@PROJECT_NAME@@**'s narrative documentation.

The TRUST documentation pipeline scans every Markdown file in this
`docs/content/` directory (the `Doxyfile`'s `INPUT` includes it, with
`RECURSIVE = YES`), so this is where your hand-written pages belong:
tutorials, theory notes, a developer guide, etc.

To grow the documentation:

1. Add a new `.md` file under `docs/content/`, opening it with a unique
   page id, e.g.:

   ```
   @page @@PROJECT_NAME@@_theory Theory
   ```

2. Link it into the navigation tree from `docs/content/main_page.md` (or
   from another page) with `@subpage`:

   ```
   - @subpage @@PROJECT_NAME@@_theory
   ```

To embed images, drop them in a directory of your choice and add it to
`IMAGE_PATH` in the `Doxyfile`'s user-override block (see the commented
example there), then reference them by their bare file name. Replace
this page with your own content once you start writing.
