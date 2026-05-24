"""Focused unit tests for trustify.cache — Phase 2 of the reshaping."""

import tempfile
import types
import unittest
from pathlib import Path


def _make_fake_content(entries):
    """Build a minimal stand-in for `TRAD2Content` that hashes deterministically.

    `cache.compute_content_hash` is documented to hash the canonical TRAD2
    text serialization, so a stub with a `toTRAD2Text()` method is enough
    for these tests — we don't need a real schema.
    """
    obj = types.SimpleNamespace()
    obj.toTRAD2Text = lambda: "\n".join(entries) + "\n"
    return obj


class TestContentHash(unittest.TestCase):
    def test_hash_is_deterministic(self):
        from trustify.cache import compute_content_hash

        c1 = _make_fake_content(["a", "b", "c"])
        c2 = _make_fake_content(["a", "b", "c"])
        self.assertEqual(
            compute_content_hash(c1, "1.2.3"),
            compute_content_hash(c2, "1.2.3"),
        )

    def test_hash_changes_when_content_changes(self):
        from trustify.cache import compute_content_hash

        c1 = _make_fake_content(["a", "b"])
        c2 = _make_fake_content(["a", "b", "c"])
        self.assertNotEqual(
            compute_content_hash(c1, "1.2.3"),
            compute_content_hash(c2, "1.2.3"),
        )

    def test_hash_changes_when_version_changes(self):
        from trustify.cache import compute_content_hash

        c = _make_fake_content(["a"])
        self.assertNotEqual(
            compute_content_hash(c, "1.2.3"),
            compute_content_hash(c, "1.2.4"),
        )

    def test_hash_is_short_hex_string(self):
        from trustify.cache import compute_content_hash

        h = compute_content_hash(_make_fake_content(["x"]), "v")
        self.assertEqual(len(h), 16)
        int(h, 16)  # raises ValueError if not hex


class TestCacheDirResolution(unittest.TestCase):
    def test_cache_dir_uses_default_root_when_no_override(self):
        from trustify.cache import cache_dir_for

        with tempfile.TemporaryDirectory() as tmp_home:
            d = cache_dir_for(
                _make_fake_content(["a"]), "1.0.0", override=None, root=Path(tmp_home) / ".cache" / "trustify"
            )
            self.assertTrue(d.exists())
            self.assertTrue(str(d).startswith(str(tmp_home)))

    def test_cache_dir_uses_override(self):
        from trustify.cache import cache_dir_for

        with tempfile.TemporaryDirectory() as tmp:
            override = Path(tmp) / "custom_out"
            d = cache_dir_for(_make_fake_content(["a"]), "1.0.0", override=override)
            self.assertEqual(d, override)
            self.assertTrue(d.exists())

    def test_is_cached_detects_both_files_present(self):
        from trustify.cache import is_cached

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            self.assertFalse(is_cached(d))
            # Pyd module name is digest-suffixed; the bare
            # `trustify_gen_pyd.py` from old caches no longer counts.
            (d / "trustify_gen_pyd_deadbeefcafebabe.py").write_text("# stub\n")
            self.assertFalse(is_cached(d))
            (d / "trustify_gen.py").write_text("# stub\n")
            # provenance.json is the completion sentinel — generate_schema
            # writes it last, after every other file is fully on disk.
            # Without it the entry is considered an interrupted write.
            self.assertFalse(is_cached(d))
            (d / "provenance.json").write_text("{}\n")
            self.assertTrue(is_cached(d))

    def test_is_cached_false_when_provenance_missing(self):
        """A partial write (gen + pyd files but no provenance.json) must
        NOT be treated as a cache hit. Reproduces audit finding 2.1:
        without the sentinel, an interrupted generate_schema's leftovers
        looked cached and the next run silently imported half-written
        files."""
        from trustify.cache import is_cached

        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "trustify_gen.py").write_text("# stub\n")
            (d / "trustify_gen_pyd_deadbeefcafebabe.py").write_text("# stub\n")
            self.assertFalse(is_cached(d))


class TestAtomicWriteDir(unittest.TestCase):
    """Helper used by generate_schema to write a cache entry without
    ever exposing a half-written final directory. Closes audit
    finding 2.1: a SIGKILL between the schema's first and last
    write must not leave the final cache path containing partial
    files. The contract:

      - The caller yields work into a sibling temp dir.
      - On normal exit, the temp dir is atomically renamed onto the
        final path (`os.replace` on POSIX = a single inode swap).
      - On exception, the temp dir is removed and the final path is
        left untouched.
      - During the `with` block, the final path does NOT yet exist
        (so concurrent readers either see the old entry or the
        complete new entry — never a torn write).
    """

    def test_atomic_write_dir_renames_temp_to_final_on_success(self):
        from trustify.cache import atomic_write_dir

        with tempfile.TemporaryDirectory() as tmp:
            final = Path(tmp) / "abc123"
            with atomic_write_dir(final) as out:
                # Final must not exist during the write window — the
                # whole point of "atomic" is that readers either see
                # nothing or see the complete entry.
                self.assertFalse(final.exists())
                (out / "trustify_gen.py").write_text("# real content\n")
            self.assertTrue(final.is_dir())
            self.assertEqual((final / "trustify_gen.py").read_text(), "# real content\n")

    def test_atomic_write_dir_cleans_up_temp_on_exception(self):
        """If the writer raises mid-write, the final dir must not be
        created AND no `.tmp.<pid>` sibling must survive."""
        from trustify.cache import atomic_write_dir

        with tempfile.TemporaryDirectory() as tmp:
            final = Path(tmp) / "abc123"
            with self.assertRaises(RuntimeError), atomic_write_dir(final) as out:
                (out / "partial.py").write_text("# half-written\n")
                raise RuntimeError("simulated crash")
            self.assertFalse(final.exists())
            # No `.abc123.tmp.<pid>` sibling left behind either.
            leftovers = [p for p in Path(tmp).iterdir() if p.name.startswith(".abc123.tmp.")]
            self.assertEqual(leftovers, [])

    def test_atomic_write_dir_replaces_preexisting_partial_dir(self):
        """If a previous crash left a partial entry at `final_dir`, the
        next atomic write replaces it cleanly. The is_cached gate
        upstream guarantees that any pre-existing dir we land on
        is_NOT a complete entry; we still must not refuse to write."""
        from trustify.cache import atomic_write_dir

        with tempfile.TemporaryDirectory() as tmp:
            final = Path(tmp) / "abc123"
            final.mkdir()
            (final / "trustify_gen.py").write_text("# leftover from crash\n")
            with atomic_write_dir(final) as out:
                (out / "trustify_gen.py").write_text("# fresh\n")
            self.assertEqual((final / "trustify_gen.py").read_text(), "# fresh\n")


class TestCacheLock(unittest.TestCase):
    """Cross-process advisory lock on the cache root. SHARED for the
    read path (LSP / is_cached + import_parser_module), EXCLUSIVE for
    the write path (generate_schema + trim). Closes audit finding 2.2:
    one process's trim must never `rmtree` the directory another
    process is currently importing from."""

    def test_cache_lock_shared_creates_and_keeps_lockfile(self):
        from trustify.cache import cache_lock

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lockfile = root / ".lock"
            self.assertFalse(lockfile.exists())
            with cache_lock("shared", root=root):
                # The lockfile is created lazily on first use.
                self.assertTrue(lockfile.is_file())
            # We deliberately leave the lockfile behind — removing it
            # would race with another process trying to acquire.
            self.assertTrue(lockfile.is_file())

    @staticmethod
    def _spawn_lock_holder(root: Path, mode: str, hold_seconds: float, ready_marker: Path):
        """Spawn a child that acquires `mode` lock, touches `ready_marker`,
        sleeps `hold_seconds`, then releases. Returns the Popen handle."""
        import subprocess
        import sys

        code = (
            "import time\n"
            "from pathlib import Path\n"
            "from trustify.cache import cache_lock\n"
            f'with cache_lock("{mode}", root=Path({str(root)!r})):\n'
            f"    Path({str(ready_marker)!r}).write_text('')\n"
            f"    time.sleep({hold_seconds!r})\n"
        )
        return subprocess.Popen([sys.executable, "-c", code])

    @staticmethod
    def _wait_for_marker(marker: Path, timeout: float = 5.0) -> None:
        import time

        deadline = time.time() + timeout
        while time.time() < deadline:
            if marker.exists():
                return
            time.sleep(0.02)
        raise AssertionError(f"child process never created marker {marker} within {timeout}s")

    def test_exclusive_lock_blocks_concurrent_exclusive(self):
        """While process A holds EXCLUSIVE, B's EXCLUSIVE attempt blocks
        until A releases. Reproduces the auto-trim-vs-write race from
        audit 2.2: two concurrent generate_schema calls cannot both be
        in the write+trim region simultaneously."""
        import time

        from trustify.cache import cache_lock

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ready = root / "READY"
            child = self._spawn_lock_holder(root, "exclusive", hold_seconds=0.8, ready_marker=ready)
            try:
                self._wait_for_marker(ready)
                # Child holds the lock. Parent's attempt must block
                # until the child releases ~0.8s after it acquired.
                start = time.time()
                with cache_lock("exclusive", root=root):
                    elapsed = time.time() - start
                self.assertGreater(elapsed, 0.2, "exclusive cache lock did not block on concurrent exclusive")
            finally:
                child.wait()

    def test_exclusive_lock_blocks_concurrent_shared(self):
        """An EXCLUSIVE writer blocks SHARED readers — the LSP can't
        import a cache entry while generate_schema is mid-write+trim."""
        import time

        from trustify.cache import cache_lock

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ready = root / "READY"
            child = self._spawn_lock_holder(root, "exclusive", hold_seconds=0.8, ready_marker=ready)
            try:
                self._wait_for_marker(ready)
                start = time.time()
                with cache_lock("shared", root=root):
                    elapsed = time.time() - start
                self.assertGreater(elapsed, 0.2, "shared cache lock did not wait for exclusive writer")
            finally:
                child.wait()

    def test_shared_holder_blocks_concurrent_exclusive(self):
        """The LSP read-path scenario: while an LSP process holds SHARED
        across is_cached + import_parser_module, a concurrent
        generate_schema call (EXCLUSIVE for its write + trim) must
        wait. Without this, the LSP's just-imported cache directory
        could be `rmtree`'d under it by an auto-trim."""
        import time

        from trustify.cache import cache_lock

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ready = root / "READY"
            child = self._spawn_lock_holder(root, "shared", hold_seconds=0.8, ready_marker=ready)
            try:
                self._wait_for_marker(ready)
                start = time.time()
                with cache_lock("exclusive", root=root):
                    elapsed = time.time() - start
                self.assertGreater(elapsed, 0.2, "exclusive must wait for active shared holders")
            finally:
                child.wait()

    def test_shared_locks_coexist(self):
        """Multiple LSPs reading the cache simultaneously must NOT
        serialise — the read path is the hot path."""
        import time

        from trustify.cache import cache_lock

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ready = root / "READY"
            child = self._spawn_lock_holder(root, "shared", hold_seconds=0.8, ready_marker=ready)
            try:
                self._wait_for_marker(ready)
                start = time.time()
                with cache_lock("shared", root=root):
                    elapsed = time.time() - start
                # Should acquire almost instantly — shared+shared compat.
                self.assertLess(elapsed, 0.2, "two shared cache locks must coexist without blocking")
            finally:
                child.wait()


class TestCacheRootEnv(unittest.TestCase):
    """The cache root honours $TRUSTIFY_CACHE_DIR.

    Lets HPC sites with a quota'd / read-only $HOME redirect the
    cache without code changes — and lets multi-user shared TRUST
    installs point every user at one shared cache.
    """

    @staticmethod
    def _set_env(name, value):
        """Save / overwrite / yield-back-a-restore-callable for one env var."""
        import os

        previous = os.environ.get(name)

        def restore():
            if previous is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = previous

        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value
        return restore

    def test_cache_dir_honours_trustify_cache_dir_env(self):
        from trustify.cache import cache_dir_for

        with tempfile.TemporaryDirectory() as tmp:
            restore = self._set_env("TRUSTIFY_CACHE_DIR", tmp)
            try:
                d = cache_dir_for(_make_fake_content(["a"]), "1.0.0")
                self.assertEqual(Path(d).parent, Path(tmp))
                self.assertTrue(d.exists())
            finally:
                restore()

    def test_default_cache_root_falls_back_when_env_unset(self):
        from trustify.cache import default_cache_root

        restore = self._set_env("TRUSTIFY_CACHE_DIR", None)
        try:
            root = default_cache_root()
        finally:
            restore()
        # Fallback is `$HOME/.cache/trustify`; we only assert the suffix
        # to stay portable across machines / CI homes.
        self.assertEqual(root.parts[-2:], (".cache", "trustify"))

    def test_default_cache_root_falls_back_when_env_is_blank(self):
        # Whitespace-only env value is treated as "unset" rather than
        # silently landing the cache in the CWD.
        from trustify.cache import default_cache_root

        restore = self._set_env("TRUSTIFY_CACHE_DIR", "   ")
        try:
            root = default_cache_root()
        finally:
            restore()
        self.assertEqual(root.parts[-2:], (".cache", "trustify"))

    def test_default_cache_root_expands_tilde(self):
        # HPC users write `TRUSTIFY_CACHE_DIR=~/scratch/...` in env
        # scripts; the leading `~` must expand.
        import os

        from trustify.cache import default_cache_root

        restore = self._set_env("TRUSTIFY_CACHE_DIR", "~/some_trustify_cache")
        try:
            root = default_cache_root()
        finally:
            restore()
        self.assertEqual(root, Path(os.path.expanduser("~/some_trustify_cache")))

    def test_iter_entries_honours_env(self):
        # `iter_entries(root=None)` is the entrypoint the CLI uses for
        # `trustify cache list`; it must follow the env var too.
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            _make_entry(tmp_path, "deadbeefcafebabe", provenance={"projects": {"x": "/y"}})
            restore = self._set_env("TRUSTIFY_CACHE_DIR", tmp)
            try:
                entries = iter_entries()
            finally:
                restore()
            self.assertEqual([e.id for e in entries], ["deadbeefcafebabe"])


def _make_entry(root: Path, entry_id: str, *, provenance: dict | None = None, extra_files: dict | None = None) -> Path:
    """Build a fake cache entry directory `root/entry_id`.

    Always writes the three is_cached sentinels:
      - `trustify_gen.py`,
      - `trustify_gen_pyd_<entry_id>.py` (production code names the
        pyd file with the same content digest as the entry directory),
      - `provenance.json` (with `provenance` content when supplied,
        else `{}` — empty content keeps `_load_entry` returning None
        fields for the "legacy / unknown provenance" test cases).
    `extra_files` are written verbatim into the entry root.
    """
    import json

    d = root / entry_id
    d.mkdir(parents=True, exist_ok=True)
    (d / "trustify_gen.py").write_text("# stub gen\n")
    (d / f"trustify_gen_pyd_{entry_id}.py").write_text("# stub pyd\n")
    (d / "provenance.json").write_text(
        json.dumps(provenance if provenance is not None else {}, indent=2, sort_keys=True) + "\n"
    )
    for name, content in (extra_files or {}).items():
        target = d / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content)
    return d


class TestIterEntries(unittest.TestCase):
    def test_empty_when_root_missing(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "nope"
            self.assertEqual(iter_entries(missing), [])

    def test_empty_when_no_cache_dirs(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            # A non-cache directory next to nothing.
            (Path(tmp) / "notacache").mkdir()
            self.assertEqual(iter_entries(Path(tmp)), [])

    def test_skips_dirs_without_generated_files(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Partial: only one of the two required files.
            (root / "partial").mkdir()
            (root / "partial" / "trustify_gen.py").write_text("x")
            _make_entry(root, "fullentry01")
            entries = iter_entries(root)
            self.assertEqual([e.id for e in entries], ["fullentry01"])

    def test_reads_provenance_fields(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "ab12cd34",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1.2.3",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            (e,) = iter_entries(root)
            self.assertEqual(e.id, "ab12cd34")
            self.assertEqual(e.created_at, "2026-05-20T18:00:00Z")
            self.assertEqual(e.trustify_version, "1.2.3")
            self.assertEqual(e.projects, {"trust": "/abs/trust"})
            self.assertGreater(e.size_bytes, 0)

    def test_entry_without_provenance_shows_none(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "ab12cd34")  # no provenance
            (e,) = iter_entries(root)
            self.assertIsNone(e.created_at)
            self.assertIsNone(e.trustify_version)
            self.assertEqual(e.projects, {})
            self.assertIsNone(e.generation_time_seconds)
            self.assertEqual(e.generation_steps_seconds, {})

    def test_reads_generation_timings(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "ab12cd34",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1.2.3",
                    "projects": {},
                    "generation_time_seconds": 12.34,
                    "generation_steps_seconds": {
                        "scan_sources": 8.0,
                        "write_trad2": 0.1,
                        "generate_pyd_and_pars": 4.24,
                    },
                },
            )
            (e,) = iter_entries(root)
            self.assertEqual(e.generation_time_seconds, 12.34)
            self.assertEqual(
                e.generation_steps_seconds,
                {
                    "scan_sources": 8.0,
                    "write_trad2": 0.1,
                    "generate_pyd_and_pars": 4.24,
                },
            )

    def test_malformed_timings_are_ignored(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "ab12cd34",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {},
                    # Non-numeric — must be dropped silently.
                    "generation_time_seconds": "nope",
                    "generation_steps_seconds": {
                        "scan_sources": "fast",
                        "valid": 1.5,
                    },
                },
            )
            (e,) = iter_entries(root)
            self.assertIsNone(e.generation_time_seconds)
            self.assertEqual(e.generation_steps_seconds, {"valid": 1.5})

    def test_sorted_newest_first(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "older123",
                provenance={
                    "generated_at": "2025-01-01T00:00:00Z",
                    "trustify_version": "1.0",
                    "projects": {},
                },
            )
            _make_entry(
                root,
                "newer456",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1.2",
                    "projects": {},
                },
            )
            ids = [e.id for e in iter_entries(root)]
            self.assertEqual(ids, ["newer456", "older123"])

    def test_entries_without_timestamp_sort_last(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "missing0", provenance=None)
            _make_entry(
                root,
                "dated001",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {},
                },
            )
            ids = [e.id for e in iter_entries(root)]
            self.assertEqual(ids, ["dated001", "missing0"])

    def test_size_bytes_counts_files_recursively(self):
        from trustify.cache import iter_entries

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "ab12cd34",
                extra_files={
                    "subdir/deep.txt": "0123456789",  # 10 bytes
                },
            )
            (e,) = iter_entries(root)
            # At minimum: deep.txt (10 bytes) + the two stub files.
            self.assertGreaterEqual(e.size_bytes, 10)
            # Confirm the recursive walk descended into subdir/.
            _make_entry(Path(tmp) / "sibling", "ef99")
            (sib,) = iter_entries(Path(tmp) / "sibling")
            self.assertGreater(e.size_bytes, sib.size_bytes)


class TestResolveId(unittest.TestCase):
    def test_full_id_resolves(self):
        from trustify.cache import resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entry = _make_entry(root, "ab12cd34567890ef")
            self.assertEqual(resolve_id("ab12cd34567890ef", root=root), entry)

    def test_short_prefix_resolves(self):
        from trustify.cache import resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entry = _make_entry(root, "ab12cd34")
            self.assertEqual(resolve_id("ab12", root=root), entry)

    def test_ambiguous_prefix_raises(self):
        from trustify.cache import AmbiguousCacheEntry, resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "abc11111")
            _make_entry(root, "abc22222")
            with self.assertRaises(AmbiguousCacheEntry) as cm:
                resolve_id("abc", root=root)
            self.assertIn("abc11111", str(cm.exception))
            self.assertIn("abc22222", str(cm.exception))

    def test_not_found_raises(self):
        from trustify.cache import CacheEntryNotFound, resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "abc11111")
            with self.assertRaises(CacheEntryNotFound):
                resolve_id("ZZ", root=root)

    def test_empty_prefix_never_matches(self):
        from trustify.cache import CacheEntryNotFound, resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "abc11111")
            with self.assertRaises(CacheEntryNotFound):
                resolve_id("", root=root)

    def test_partial_dirs_excluded_from_matches(self):
        from trustify.cache import CacheEntryNotFound, resolve_id

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A dir without generated files must not match.
            (root / "abc99999").mkdir()
            with self.assertRaises(CacheEntryNotFound):
                resolve_id("abc", root=root)


class TestCopyEntry(unittest.TestCase):
    def test_copies_files(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(
                root,
                "abc11111",
                provenance={
                    "generated_at": "2026-01-01T00:00:00Z",
                    "trustify_version": "1",
                    "projects": {},
                },
            )
            dest = Path(tmp) / "dest"
            copy_entry(src, dest)
            self.assertTrue((dest / "trustify_gen.py").is_file())
            self.assertTrue((dest / "trustify_gen_pyd_abc11111.py").is_file())
            self.assertTrue((dest / "provenance.json").is_file())

    def test_creates_missing_dest_dir(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(root, "abc11111")
            dest = Path(tmp) / "new" / "nested" / "dest"
            copy_entry(src, dest)
            self.assertTrue(dest.is_dir())

    def test_filters_pycache_at_any_depth(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(
                root,
                "abc11111",
                extra_files={
                    "__pycache__/foo.cpython-312.pyc": "bytecode",
                    "nested/__pycache__/bar.cpython-312.pyc": "bytecode",
                    "nested/stay.txt": "kept",
                },
            )
            dest = Path(tmp) / "dest"
            copy_entry(src, dest)
            self.assertFalse((dest / "__pycache__").exists())
            self.assertFalse((dest / "nested" / "__pycache__").exists())
            self.assertTrue((dest / "nested" / "stay.txt").is_file())

    def test_filters_pyc_at_top_level(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(
                root,
                "abc11111",
                extra_files={
                    "stray.pyc": "bytecode",
                },
            )
            dest = Path(tmp) / "dest"
            copy_entry(src, dest)
            self.assertFalse((dest / "stray.pyc").exists())

    def test_non_empty_dest_raises(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(root, "abc11111")
            dest = Path(tmp) / "dest"
            dest.mkdir()
            (dest / "blocker").write_text("x")
            with self.assertRaises(FileExistsError):
                copy_entry(src, dest)

    def test_dest_is_a_file_raises(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(root, "abc11111")
            dest = Path(tmp) / "afile"
            dest.write_text("not a dir")
            with self.assertRaises(FileExistsError):
                copy_entry(src, dest)

    def test_empty_dest_ok(self):
        from trustify.cache import copy_entry

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = _make_entry(root, "abc11111")
            dest = Path(tmp) / "dest"
            dest.mkdir()  # exists, empty
            copy_entry(src, dest)
            self.assertTrue((dest / "trustify_gen.py").is_file())


class TestTrimSameProjectSet(unittest.TestCase):
    """`trim_same_project_set` is the auto-trim hook fired by
    api.generate_schema after writing a new entry."""

    def test_evicts_other_entries_with_same_projects(self):
        from trustify.cache import trim_same_project_set

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keeper = _make_entry(
                root,
                "kept11111",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "stale111",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            deleted = trim_same_project_set(keeper, root=root)
            self.assertEqual(deleted, ["stale111"])
            self.assertTrue(keeper.exists())
            self.assertFalse((root / "stale111").exists())

    def test_preserves_entries_with_different_projects(self):
        from trustify.cache import trim_same_project_set

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keeper = _make_entry(
                root,
                "kept11111",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "baltik11",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust", "baltik_a": "/abs/b"},
                },
            )
            self.assertEqual(trim_same_project_set(keeper, root=root), [])
            self.assertTrue((root / "baltik11").exists())

    def test_preserves_entries_with_different_paths(self):
        # Same project names but different absolute paths (multi-user
        # on shared $HOME / two clones in different dirs) → keep both.
        from trustify.cache import trim_same_project_set

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keeper = _make_entry(
                root,
                "kept11111",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "other111",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/different/path"},
                },
            )
            self.assertEqual(trim_same_project_set(keeper, root=root), [])
            self.assertTrue((root / "other111").exists())

    def test_no_match_when_keep_has_no_provenance(self):
        from trustify.cache import trim_same_project_set

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keeper = _make_entry(root, "kept11111")
            _make_entry(
                root,
                "stale111",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            # The keeper has no projects → no group to dedup.
            self.assertEqual(trim_same_project_set(keeper, root=root), [])
            self.assertTrue((root / "stale111").exists())

    def test_skips_entries_without_provenance(self):
        from trustify.cache import trim_same_project_set

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keeper = _make_entry(
                root,
                "kept11111",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(root, "noprov11")  # pre-refactor entry, no provenance
            self.assertEqual(trim_same_project_set(keeper, root=root), [])
            self.assertTrue((root / "noprov11").exists())


class TestTrimAllDuplicates(unittest.TestCase):
    """`trim_all_duplicates` powers the `trustify cache trim` subcommand."""

    def test_keeps_newest_per_group(self):
        from trustify.cache import trim_all_duplicates

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "trust_a1",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "trust_a2",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "balt_b1",
                provenance={
                    "generated_at": "2026-05-20T17:00:00Z",
                    "trustify_version": "1",
                    "projects": {"baltik_a": "/abs/b"},
                },
            )
            _make_entry(
                root,
                "balt_b2",
                provenance={
                    "generated_at": "2026-05-18T17:00:00Z",
                    "trustify_version": "1",
                    "projects": {"baltik_a": "/abs/b"},
                },
            )
            result = trim_all_duplicates(root=root)
            # Two groups, each with one keeper + one eviction.
            self.assertEqual(set(result.keys()), {"trust_a1", "balt_b1"})
            self.assertEqual(result["trust_a1"], ["trust_a2"])
            self.assertEqual(result["balt_b1"], ["balt_b2"])
            self.assertFalse((root / "trust_a2").exists())
            self.assertFalse((root / "balt_b2").exists())
            self.assertTrue((root / "trust_a1").exists())
            self.assertTrue((root / "balt_b1").exists())

    def test_no_duplicates_returns_empty(self):
        from trustify.cache import trim_all_duplicates

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(
                root,
                "trust_a1",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_entry(
                root,
                "balt_b1",
                provenance={
                    "generated_at": "2026-05-20T17:00:00Z",
                    "trustify_version": "1",
                    "projects": {"baltik_a": "/abs/b"},
                },
            )
            self.assertEqual(trim_all_duplicates(root=root), {})

    def test_entries_without_provenance_are_ignored(self):
        from trustify.cache import trim_all_duplicates

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_entry(root, "noprov11")
            _make_entry(root, "noprov22")
            # Two no-provenance entries are never grouped together
            # (we don't have data to compare them on).
            self.assertEqual(trim_all_duplicates(root=root), {})
            self.assertTrue((root / "noprov11").exists())
            self.assertTrue((root / "noprov22").exists())


class TestShouldAutoTrim(unittest.TestCase):
    def test_default_is_true(self):
        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            # No .trustify.json anywhere walking up from tmp → default True.
            # We deliberately use a temp tree to dodge any real config
            # that might live above this checkout.
            self.assertTrue(should_auto_trim(start_dir=Path(tmp)))

    def test_env_var_disables(self):
        import os

        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            old = os.environ.get("TRUSTIFY_NO_AUTO_TRIM")
            try:
                os.environ["TRUSTIFY_NO_AUTO_TRIM"] = "1"
                self.assertFalse(should_auto_trim(start_dir=Path(tmp)))
            finally:
                if old is None:
                    os.environ.pop("TRUSTIFY_NO_AUTO_TRIM", None)
                else:
                    os.environ["TRUSTIFY_NO_AUTO_TRIM"] = old

    def test_env_var_any_non_empty_value_disables(self):
        """Audit 1.6: `TRUSTIFY_NO_AUTO_TRIM` follows the conventional
        env-var-presence semantic (same as `NO_COLOR` elsewhere in
        trustify): any non-empty value disables. The old deny-list
        for `"0"` / `"false"` / `"no"` was a surprise — users
        reasonably expect `TRUSTIFY_NO_AUTO_TRIM=0` to "set" the
        variable and therefore disable, not the reverse.
        """
        import os

        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            for val in ("1", "0", "false", "FALSE", "no", "yes", "anything", "   "):
                old = os.environ.get("TRUSTIFY_NO_AUTO_TRIM")
                try:
                    os.environ["TRUSTIFY_NO_AUTO_TRIM"] = val
                    self.assertFalse(should_auto_trim(start_dir=Path(tmp)), f"value {val!r} should disable")
                finally:
                    if old is None:
                        os.environ.pop("TRUSTIFY_NO_AUTO_TRIM", None)
                    else:
                        os.environ["TRUSTIFY_NO_AUTO_TRIM"] = old

    def test_env_var_empty_string_treated_as_unset(self):
        """`TRUSTIFY_NO_AUTO_TRIM=""` is falsy (matching NO_COLOR's
        behaviour in `_color.py`): the variable is "set to nothing",
        which doesn't carry the user's intent to disable. Equivalent
        to unset → default (enabled).
        """
        import os

        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            old = os.environ.get("TRUSTIFY_NO_AUTO_TRIM")
            try:
                os.environ["TRUSTIFY_NO_AUTO_TRIM"] = ""
                self.assertTrue(should_auto_trim(start_dir=Path(tmp)))
            finally:
                if old is None:
                    os.environ.pop("TRUSTIFY_NO_AUTO_TRIM", None)
                else:
                    os.environ["TRUSTIFY_NO_AUTO_TRIM"] = old

    def test_config_file_disables(self):
        import json as _json

        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            cfg = Path(tmp) / ".trustify.json"
            cfg.write_text(_json.dumps({"cache": {"auto_trim": False}}))
            self.assertFalse(should_auto_trim(start_dir=Path(tmp)))

    def test_config_file_explicit_true_keeps_default(self):
        import json as _json

        from trustify.cache import should_auto_trim

        with tempfile.TemporaryDirectory() as tmp:
            cfg = Path(tmp) / ".trustify.json"
            cfg.write_text(_json.dumps({"cache": {"auto_trim": True}}))
            self.assertTrue(should_auto_trim(start_dir=Path(tmp)))


class TestEntrySizeHandlesBrokenSymlinks(unittest.TestCase):
    """Audit 5.5: `_entry_size` walks with `followlinks=False`, but
    inside the loop it used `os.path.getsize` — which DOES follow
    symlinks. A broken file symlink in the cache dir then raised
    `OSError`, the OSError-tolerant `continue` swallowed it, and
    `trustify cache list` understated the entry size by the symlink
    slot's own bytes. Switching to `os.lstat(...).st_size` keeps the
    semantic consistent with the walk (count the link slot itself,
    never the target), so a broken symlink contributes its actual
    on-disk footprint rather than zero.
    """

    def test_broken_symlink_contributes_its_own_size(self):
        import os

        from trustify.cache import _entry_size

        with tempfile.TemporaryDirectory() as tmp:
            entry = Path(tmp)
            (entry / "regular.txt").write_bytes(b"x" * 100)
            broken_target = "/nonexistent_target_with_a_known_length"
            os.symlink(broken_target, entry / "dangling")
            total = _entry_size(entry)
            # Regular file (100 B) + symlink slot (len of target path)
            self.assertEqual(total, 100 + len(broken_target))

    def test_unbroken_symlink_contributes_its_own_size_not_target(self):
        # Consistency check: even a live symlink reports its own bytes,
        # not the target's — same semantic as the walk's followlinks=False.
        import os

        from trustify.cache import _entry_size

        with tempfile.TemporaryDirectory() as tmp:
            entry = Path(tmp)
            target = entry / "target.txt"
            target.write_bytes(b"y" * 500)
            link = entry / "alias"
            os.symlink(str(target), link)
            total = _entry_size(entry)
            # 500 (target file) + len(target absolute path) (symlink slot).
            self.assertEqual(total, 500 + len(str(target)))


if __name__ == "__main__":
    unittest.main()
