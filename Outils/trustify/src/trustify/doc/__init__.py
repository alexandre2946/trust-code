"""Markdown keyword-reference generation for `trustify generate_markdown`.

Public surface:
- ``trustify.doc.extras``    — ``ExtrasIndex`` (per-project extras scan
                               + lookup + usage tracking).
- ``trustify.doc.origin``    — ``OriginClassifier`` (per-class origin
                               classifier, --projects-driven).
- ``trustify.doc.render``    — ``render_keyword_card`` (one Markdown
                               card per keyword class).
- ``trustify.doc.generator`` — ``DocGenerator`` (orchestrator: imports
                               the generated pydantic module, groups
                               classes by (origin, family), emits
                               hub/branch/index pages).

End-user entry point is ``trustify.api.generate_markdown``; the CLI hits it via
``trustify generate_markdown`` (see ``cli/__init__.py``).
"""

from trustify.doc.extras import ExtrasIndex
from trustify.doc.generator import DocGenerator
from trustify.doc.origin import OriginClassifier

__all__ = ["DocGenerator", "ExtrasIndex", "OriginClassifier"]
