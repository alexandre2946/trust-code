#!/usr/bin/env python3
"""Tests for refactor_ptr.refactor().

We test two modes:

  * NO-WRAP mode (`wrap_bool=False`): the rewrite always emits a bare
    receiver for `non_nul()`. Correct when the pointer type has an
    implicit `operator bool`.

  * WRAP mode (`wrap_bool=True`, default): the rewrite wraps the bare
    receiver in `bool(...)` unless the context is provably boolean
    (operand of `!`, `&&`, `||`, ternary `?`, or inside the condition
    parens of `if`/`while`/`for`). Correct for `explicit operator bool`.
"""
import sys
sys.path.insert(0, '/home/claude/refactor')
from refactor_ptr import refactor


# ---------------------------------------------------------------------------
# NO-WRAP mode tests — the original behavior.
# ---------------------------------------------------------------------------
NO_WRAP_CASES = [
    # basic
    ("ptr.est_nul()",          "!ptr"),
    ("ptr.non_nul()",          "ptr"),

    # leading ! (collapse/absorb)
    ("!ptr.est_nul()",         "ptr"),
    ("!ptr.non_nul()",         "!ptr"),
    ("! ptr.est_nul()",        "ptr"),
    ("!  ptr.non_nul()",       "!ptr"),
    ("!!ptr.est_nul()",        "!!!ptr"),

    # member access and arrow
    ("obj.ptr.est_nul()",              "!obj.ptr"),
    ("obj->ptr.est_nul()",             "!obj->ptr"),
    ("obj->member->ptr.non_nul()",     "obj->member->ptr"),
    ("a.b.c.d.est_nul()",              "!a.b.c.d"),
    ("Foo::Bar::instance.non_nul()",   "Foo::Bar::instance"),

    # function calls in the receiver chain
    ("getPtr().est_nul()",             "!getPtr()"),
    ("getPtr(a, b).non_nul()",         "getPtr(a, b)"),
    ("obj.getChild(1,2).est_nul()",    "!obj.getChild(1,2)"),
    ("foo(bar(baz)).non_nul()",        "foo(bar(baz))"),

    # subscripts
    ("arr[i].est_nul()",               "!arr[i]"),
    ("arr[i+1][j].non_nul()",          "arr[i+1][j]"),
    ("map[\"key\"].est_nul()",         "!map[\"key\"]"),

    # mixed chain
    ("obj->items[i]->get().est_nul()", "!obj->items[i]->get()"),
    ("!obj->items[i]->get().est_nul()", "obj->items[i]->get()"),

    # boolean expressions
    ("if (x.est_nul()) { }",           "if (!x) { }"),
    ("if (x.non_nul()) { }",           "if (x) { }"),
    ("if (x.est_nul() && y.non_nul()) { }", "if (!x && y) { }"),
    ("if (!x.est_nul() || y.non_nul()) { }", "if (x || y) { }"),
    ("return x.est_nul();",            "return !x;"),
    ("return !x.non_nul();",           "return !x;"),

    # ternary
    ("a = p.non_nul() ? p->val : 0;",  "a = p ? p->val : 0;"),
    ("a = !p.est_nul() ? p->val : 0;", "a = p ? p->val : 0;"),

    # whitespace variants
    ("ptr . est_nul ( )",              "!ptr"),
    ("ptr.est_nul(  )",                "!ptr"),

    # strings/comments untouched
    ('const char* s = "ptr.est_nul()";', 'const char* s = "ptr.est_nul()";'),
    ("// ptr.est_nul() stays\nreal.est_nul();",
     "// ptr.est_nul() stays\n!real;"),
    ("/* ptr.non_nul() */ x.non_nul();", "/* ptr.non_nul() */ x;"),
    ('s = "a\\"b.est_nul()"; x.est_nul();',
     's = "a\\"b.est_nul()"; !x;'),
    ('auto s = R"(ptr.est_nul())"; x.non_nul();',
     'auto s = R"(ptr.est_nul())"; x;'),
    ('auto s = R"xy(p.est_nul())xy"; x.est_nul();',
     'auto s = R"xy(p.est_nul())xy"; !x;'),
    ("char c = '.'; x.est_nul();",     "char c = '.'; !x;"),

    # multiple occurrences on one line
    ("if (a.est_nul() || b.est_nul()) break;",
     "if (!a || !b) break;"),
    ("x.non_nul() && y.non_nul() && z.non_nul()",
     "x && y && z"),

    # nested calls
    ("foo(x.non_nul()).est_nul()",     "!foo(x)"),
    ("foo(x.est_nul(), y.non_nul()).non_nul()",
     "foo(!x, y)"),

    # != not mistaken for leading !
    ("a != b.est_nul()",               "a != !b"),
    ("a!=b.non_nul()",                 "a!=b"),

    # non-negation left context
    ("x = p.est_nul();",               "x = !p;"),
    ("(p.est_nul())",                  "(!p)"),
    ("!(p.est_nul())",                 "!(!p)"),

    # templates, this->
    ("obj.get<Foo>().est_nul()",       "!obj.get<Foo>()"),
    ("this->m_ptr.est_nul()",          "!this->m_ptr"),
    ("!this->m_ptr.non_nul()",         "!this->m_ptr"),

    # must not rewrite unrelated methods
    ("obj.est_nul_other()",            "obj.est_nul_other()"),

    # keywords and surrounding contexts
    ("if(p.est_nul())",                "if(!p)"),
    ("while (iter.non_nul()) {",       "while (iter) {"),

    # multi-line expression
    ("if (foo\n    .est_nul())",       "if (!foo)"),
    ("auto x = foo\n    ->bar\n    .non_nul();",
     "auto x = foo\n    ->bar;"),

    # negation in larger expression
    ("x = a && !b.est_nul();",         "x = a && b;"),
    ("x = !a || !b.non_nul();",        "x = !a || !b;"),
    ("x = 3 != y.non_nul();",          "x = 3 != y;"),

    # chain with templates and calls
    ("obj.get<Foo, Bar>().child.est_nul()",
     "!obj.get<Foo, Bar>().child"),
]


# ---------------------------------------------------------------------------
# WRAP mode tests — default behavior for `explicit operator bool`.
# ---------------------------------------------------------------------------
WRAP_CASES = [
    # ==== non_nul in NON-boolean context: must be wrapped ====

    # The motivating bug report.
    ("bool is_coupled() const { return pbc_.non_nul(); }",
     "bool is_coupled() const { return bool(pbc_); }"),

    # Assignment
    ("bool ok = ptr.non_nul();",         "bool ok = bool(ptr);"),
    ("auto x = ptr.non_nul();",          "auto x = bool(ptr);"),

    # Function argument
    ("log(ptr.non_nul());",              "log(bool(ptr));"),
    ("log(a, ptr.non_nul(), b);",        "log(a, bool(ptr), b);"),

    # Bare expression statement
    ("ptr.non_nul();",                   "bool(ptr);"),

    # != comparison: bool(...) still needed (explicit bool doesn't convert
    # implicitly through operator==/operator!=).
    ("x = 3 != y.non_nul();",            "x = 3 != bool(y);"),

    # ==== non_nul in BOOLEAN context: must NOT be wrapped ====

    ("if (ptr.non_nul()) {}",            "if (ptr) {}"),
    ("while (ptr.non_nul()) {}",         "while (ptr) {}"),
    ("for (; ptr.non_nul(); )",          "for (; ptr; )"),
    ("if (a.non_nul() && b.non_nul())",  "if (a && b)"),
    ("if (a.non_nul() || b)",            "if (a || b)"),
    ("x = a && ptr.non_nul();",          "x = a && ptr;"),
    ("x = ptr.non_nul() && a;",          "x = ptr && a;"),
    ("x = ptr.non_nul() || other;",      "x = ptr || other;"),
    ("x = ptr.non_nul() ? 1 : 0;",       "x = ptr ? 1 : 0;"),
    ("!ptr.non_nul()",                   "!ptr"),
    ("if (!ptr.non_nul())",              "if (!ptr)"),
    ("assert(ptr.non_nul());",           "assert(ptr);"),

    # ==== est_nul never wraps (the ! forces conversion) ====

    ("bool is_null() const { return p.est_nul(); }",
     "bool is_null() const { return !p; }"),
    ("bool x = ptr.est_nul();",          "bool x = !ptr;"),
    ("log(ptr.est_nul());",              "log(!ptr);"),
    ("if (ptr.est_nul())",               "if (!ptr)"),

    # ==== !ptr.est_nul() -> ptr (double-negation collapse) ====
    # The resulting bare receiver needs wrapping in non-boolean context.

    ("return !p.est_nul();",             "return bool(p);"),
    ("auto x = !p.est_nul();",           "auto x = bool(p);"),
    ("if (!p.est_nul())",                "if (p)"),
    ("while (!p.est_nul()) {}",          "while (p) {}"),
    ("x = !p.est_nul() && y;",           "x = p && y;"),

    # ==== complex receivers wrapped correctly ====

    ("return obj->ptr.non_nul();",       "return bool(obj->ptr);"),
    ("return getPtr().non_nul();",       "return bool(getPtr());"),
    ("bool b = arr[i].non_nul();",       "bool b = bool(arr[i]);"),
    ("return obj.get<Foo>().non_nul();", "return bool(obj.get<Foo>());"),

    # ==== for-loop: only the condition slot is boolean ====

    ("for (auto p = get(); p.non_nul(); p = next()) {}",
     "for (auto p = get(); p; p = next()) {}"),
    ("for (ptr.non_nul(); i < 10; ++i) {}",
     "for (bool(ptr); i < 10; ++i) {}"),

    # ==== strings/comments must still be untouched ====

    ('const char* s = "ptr.non_nul()"; return ptr.non_nul();',
     'const char* s = "ptr.non_nul()"; return bool(ptr);'),
]


def run_suite(name, cases, **kwargs):
    passed = 0
    failed = 0
    for i, (src, expected) in enumerate(cases):
        got, _ = refactor(src, **kwargs)
        if got == expected:
            passed += 1
        else:
            failed += 1
            print(f"FAIL [{name}] case {i}:")
            print(f"  input:    {src!r}")
            print(f"  expected: {expected!r}")
            print(f"  got:      {got!r}")
    print(f"[{name}] {passed}/{passed+failed} passed")
    return failed


def run():
    fails = 0
    fails += run_suite("no-wrap", NO_WRAP_CASES, wrap_bool=False)
    fails += run_suite("wrap",    WRAP_CASES,    wrap_bool=True)
    return 0 if fails == 0 else 1


if __name__ == '__main__':
    sys.exit(run())
