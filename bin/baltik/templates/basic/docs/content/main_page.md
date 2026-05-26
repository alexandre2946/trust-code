@mainpage @@PROJECT_NAME@@ Documentation

Welcome to the **@@PROJECT_NAME@@** documentation.

This documentation covers the full C++ API of **@@PROJECT_NAME@@** together
with its TRUST platform dependencies.

---

<!-- Hand-written narrative pages live under docs/content/ (scanned via the
     Doxyfile's INPUT). demo.md ships as a starting point — replace it with
     your own pages and list each one here with @subpage. -->
- @subpage @@PROJECT_NAME@@_demo

<!-- The 'KeywordsReference' page id below must match the @page
     declaration emitted by trustify (Outils/trustify/src/trustify/doc/
     generator.py — see `@page KeywordsReference Keywords Reference Manual`).
     If trustify ever renames this page, update the @subpage reference
     here AND in TRUST's docs/content/main_page.md. -->
- @subpage KeywordsReference

---

<h2>Quick start</h2>

Source the environment and build the project:

```bash
source env_@@PROJECT_NAME@@.sh
make optim
```

Then generate this documentation:

```bash
make docs
# or directly:
cd docs && make
```

The HTML output is written to `docs/build/html/index.html`.

<h2>Credits</h2>

**@@PROJECT_NAME@@** is built on top of the [TRUST platform](https://cea-trust-platform.github.io/),
an open-source thermohydraulic framework developed at CEA.

- Github: <https://github.com/cea-trust-platform/trust-code>
- Documentation: <https://cea-trust-platform.readthedocs.io/en/>

