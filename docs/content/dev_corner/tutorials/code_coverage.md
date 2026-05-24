@page DevTut_CodeCoverage Code coverage: targeted test selection

When you touch a class or method in your baltik (or in TRUST itself),
you usually don't need to run the full non-regression suite — only the
test cases that exercise the code path you changed.  The `trust -check`
flag has a *coverage* variant that selects test cases by their C++
symbol coverage, letting you re-run exactly the relevant ones in
seconds rather than hours.

## The -check coverage forms

In addition to running named test cases, `trust -check` accepts a
**function**, **class** or **class::method** name.  It looks up which
test cases exercise that symbol (via the coverage database built by
the nightly CI) and runs only those:

```bash
trust -check FUNCTION_NAME
trust -check CLASS_NAME
trust -check CLASS_NAME::METHOD_NAME
```

`-ctest` accepts the same forms and routes through the CTest harness
when you prefer that runner:

```bash
trust -ctest CLASS_NAME::METHOD_NAME
```

## Worked example: the RRK2 rational Runge-Kutta scheme

Say you've modified the `RRK2` class (rational Runge-Kutta order 2)
and you want to re-run only the test cases that touch it.  Find the
class in the Doxygen documentation first — that confirms the public
methods and gives you their exact spellings.  Then:

```bash
trust -check RRK2::RRK2          # only the cases hitting the constructor
trust -check RRK2                # every case hitting any RRK2 method
```

The first invocation is the laser-focused one: it re-runs only the
test cases whose execution actually entered `RRK2::RRK2`.  The second
is broader and useful as a smoke test of the class as a whole before
opening a pull request.

## When to use coverage-driven selection

- **Before a commit**, to confirm your fix doesn't regress any of the
  test cases that exercise the modified code.
- **During bisection** when narrowing down which test exposes a
  regression.
- **In CI** for fast feedback on pull-request branches before running
  the full nightly battery.

The coverage information is populated by the periodic CI runs that
collect gcov data on the test suite; if a recently-added test case
doesn't yet appear in the selection, the next coverage run will pick
it up.
