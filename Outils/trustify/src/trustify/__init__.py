try:
    from importlib.metadata import version as _pkg_version

    __version__ = _pkg_version("trustify")
except Exception:
    # Source checkout that hasn't been pip-installed, or a runtime where
    # importlib.metadata can't find the dist info — fall back to a sentinel
    # rather than crashing the import chain.
    __version__ = "0.0.0+unknown"

from trustify.api import (  # noqa: F401
    batch_check,
    check,
    format_dataset_file,
    generate_schema,
    init_config,
    load_dataset,
    modernize,
)
from trustify.core.base import Dataset_Parser  # noqa: F401
from trustify.core.misc_utilities import TrustifyInternalError, TrustifyParseError  # noqa: F401
from trustify.core.trust_parser import SourceRange, TRUSTParser, TRUSTStream  # noqa: F401
from trustify.formatter import format_dataset  # noqa: F401

# Public opt-in for programmatic callers that want CLI-style workspace
# resolution (see docs/project-resolution.md). The api.* entry points
# never walk the filesystem on their own — pass `discover_config(path)`
# results through explicitly when needed.
from trustify.projects import (  # noqa: F401
    ConfigError,
    TrustifyConfig,
    detect_baltik_root,
    discover_config,
    expand_dependencies,
)
