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
from igv_helpers import per_read_warning_note

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
