"""Tests for per_read_warning_note — the per-read filter banner text helper.

The bug being tested: the warning banner previously used a single hardcoded
string that always claimed "alt_count and VAF are re-aggregated from reads
passing the filter", regardless of whether recompute_vaf was True or False.
In the default locus-inclusion mode (recompute_vaf=False) alt_count is NOT
re-aggregated, so the message was factually wrong.

The fix extracted the note into per_read_warning_note(recompute_vaf) which
returns mode-appropriate text.

Run with: pytest app/tests/
"""

import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from igv_helpers import per_read_warning_note, insert_size_active_part

# The exact phrase that was always shown before the fix.
_OLD_HARDCODED_PHRASE = "alt_count and VAF are re-aggregated from reads passing the filter"


class TestPerReadWarningNote:

    def test_recompute_mode_mentions_reaggregation(self):
        """In recompute_vaf=True mode the note must mention re-aggregation."""
        note = per_read_warning_note(recompute_vaf=True)
        assert "re-aggregated" in note

    def test_locus_inclusion_mode_does_not_mention_reaggregation(self):
        """In the default locus-inclusion mode the note must NOT claim re-aggregation.

        This is the regression test for bug #2. Before the fix, the banner
        always contained the re-aggregation phrase even when recompute_vaf=False.
        """
        note = per_read_warning_note(recompute_vaf=False)
        assert "re-aggregated" not in note

    def test_locus_inclusion_mode_mentions_hidden_loci(self):
        """Locus-inclusion mode should tell the user that non-passing loci are hidden."""
        note = per_read_warning_note(recompute_vaf=False)
        assert "hidden" in note.lower()

    def test_recompute_mode_mentions_original_vaf(self):
        """Re-aggregation mode should reference original_vaf for comparison."""
        note = per_read_warning_note(recompute_vaf=True)
        assert "original_vaf" in note

    def test_buggy_hardcoded_string_would_be_wrong_for_locus_inclusion(self):
        """Confirm the old hardcoded phrase is incorrect for locus-inclusion mode.

        Documents why the fix was necessary: the old single string falsely
        claimed re-aggregation in the default (non-recompute) mode.
        """
        note = per_read_warning_note(recompute_vaf=False)
        assert _OLD_HARDCODED_PHRASE not in note, (
            "The old hardcoded re-aggregation phrase must not appear in "
            "locus-inclusion mode (recompute_vaf=False)"
        )

    def test_two_modes_return_different_text(self):
        """The two modes must produce distinct messages."""
        assert per_read_warning_note(True) != per_read_warning_note(False)


class TestInsertSizeActivePart:
    """Tests for insert_size_active_part — bug #3 regression + exclude mode.

    The bug: _active_parts in the warning banner never included insert size,
    so activating only the insert size filter produced "Per-read filters active ()".
    The fix adds insert_size_active_part() and appends its result when non-None.
    """

    _MIN, _MAX = 20, 500  # mirrors _IS_MIN, _IS_MAX in geac_explorer.py

    def test_returns_none_at_defaults(self):
        """At slider defaults no part is emitted — insert size filter is inactive."""
        assert insert_size_active_part(self._MIN, self._MAX, self._MIN, self._MAX) is None

    def test_returns_part_when_lo_raised(self):
        """Raising the lower bound makes the filter active."""
        part = insert_size_active_part(50, self._MAX, self._MIN, self._MAX)
        assert part is not None
        assert "insert size" in part

    def test_returns_part_when_hi_lowered(self):
        """Lowering the upper bound makes the filter active."""
        part = insert_size_active_part(self._MIN, 400, self._MIN, self._MAX)
        assert part is not None
        assert "insert size" in part

    def test_part_contains_bounds(self):
        """The emitted string should contain the active lo and hi values."""
        part = insert_size_active_part(50, 400, self._MIN, self._MAX)
        assert "50" in part
        assert "400" in part

    def test_include_mode_says_including(self):
        """Include mode (default) should say 'including only'."""
        part = insert_size_active_part(50, 400, self._MIN, self._MAX, is_exclude_mode=False)
        assert "including only" in part

    def test_exclude_mode_says_excluding(self):
        """Exclude mode should say 'excluding'."""
        part = insert_size_active_part(50, 400, self._MIN, self._MAX, is_exclude_mode=True)
        assert "excluding" in part

    def test_include_and_exclude_modes_differ(self):
        """The two modes must produce distinct strings."""
        inc = insert_size_active_part(50, 400, self._MIN, self._MAX, is_exclude_mode=False)
        exc = insert_size_active_part(50, 400, self._MIN, self._MAX, is_exclude_mode=True)
        assert inc != exc

    def test_buggy_code_would_omit_insert_size_from_banner(self):
        """Confirm the bug: before the fix, _active_parts never got an insert-size entry.

        Simulate the old code path (no insert_size_active_part call) and the
        new one side-by-side to document exactly what the fix adds.
        """
        # Old code: _active_parts built without any insert size logic
        old_active_parts = []  # insert size was simply never appended

        # New code: insert_size_active_part is called and appended when non-None
        new_active_parts = []
        part = insert_size_active_part(50, 400, self._MIN, self._MAX)
        if part is not None:
            new_active_parts.append(part)

        assert not any("insert size" in p for p in old_active_parts), (
            "Old code should have no insert-size entry (documents the bug)"
        )
        assert any("insert size" in p for p in new_active_parts), (
            "New code must include insert-size entry when filter is active"
        )
