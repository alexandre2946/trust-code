# Golden-file Formatting Corpus

Pairs of `inputs/<name>.data` + `expected/<name>.data` files that define the expected formatter output for `.data` inputs.

Consumed by the Python formatter tests (`tests/test_corpus.py` in the `trustify` package).

All cases assume the locked-in formatter parameters: 4-space indentation, no startContext.

Add a new pair here whenever you fix a formatter bug or change behavior. The corpus is the behavioral contract — everything else (unit tests, LSP integration tests) must agree with it.

## Notes on specific cases

- `template-braced-substitution.data` / `template-bare-substitution.data` — locks in that the formatter passes through Python `string.Template` placeholders (`${key}` and `$key`) unchanged. These appear in real TRUST validation datasets that get materialised at notebook-run time by `Template(...).substitute(dict)`. Excerpted from `Validation/Rapports_automatiques/Verification/Verification_codage/Diffusion_DG/src/jdd_test.data` and `.../Multiphase/canal_axi/src/jdd_pbmulti.data` respectively. The parser and LSP do not handle these (the templates are not valid TRUST syntax until substituted) but the formatter must, so users can pre-format the template before substitution.
